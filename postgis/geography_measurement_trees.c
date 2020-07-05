/**********************************************************************
 *
 * PostGIS - Spatial Types for PostgreSQL
 * http://postgis.net
 *
 * PostGIS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * PostGIS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PostGIS.  If not, see <http://www.gnu.org/licenses/>.
 *
 **********************************************************************
 *
 * ^copyright^
 *
 **********************************************************************/

#include "geography_measurement_trees.h"


/*
* Specific tree types include all the generic slots and
* their own slots for their trees. We put the implementation
* for the CircTreeGeomCache here because we can't shove
* the PgSQL specific bits of the code (fcinfo) back into
* liblwgeom, where most of the circtree logic lives.
*/
typedef struct {
	GeomCache    gcache;
	CIRC_NODE*   index;
} CircTreeGeomCache;



/**
* Builder, freeer and public accessor for cached CIRC_NODE trees
*/
static int
CircTreeBuilder(const LWGEOM* lwgeom, GeomCache* cache)
{
	CircTreeGeomCache* circ_cache = (CircTreeGeomCache*)cache;
	CIRC_NODE* tree = lwgeom_calculate_circ_tree(lwgeom);

	if ( circ_cache->index )
	{
		circ_tree_free(circ_cache->index);
		circ_cache->index = 0;
	}
	if ( ! tree )
		return LW_FAILURE;

	circ_cache->index = tree;
	return LW_SUCCESS;
}

static int
CircTreeFreer(GeomCache* cache)
{
	CircTreeGeomCache* circ_cache = (CircTreeGeomCache*)cache;
	if ( circ_cache->index )
	{
		circ_tree_free(circ_cache->index);
		circ_cache->index = 0;
		circ_cache->gcache.argnum = 0;
	}
	return LW_SUCCESS;
}

static GeomCache*
CircTreeAllocator(void)
{
	CircTreeGeomCache* cache = palloc(sizeof(CircTreeGeomCache));
	memset(cache, 0, sizeof(CircTreeGeomCache));
	return (GeomCache*)cache;
}

static GeomCacheMethods CircTreeCacheMethods =
{
	CIRC_CACHE_ENTRY,
	CircTreeBuilder,
	CircTreeFreer,
	CircTreeAllocator
};

static CircTreeGeomCache *
GetCircTreeGeomCache(FunctionCallInfo fcinfo, SHARED_GSERIALIZED *g1, SHARED_GSERIALIZED *g2)
{
	return (CircTreeGeomCache*)GetGeomCache(fcinfo, &CircTreeCacheMethods, g1, g2);
}

static int
CircTreePIP(const CIRC_NODE* tree1, const GSERIALIZED* g1, const POINT4D* in_point)
{
	int tree1_type = gserialized_get_type(g1);
	GBOX gbox1;
	GEOGRAPHIC_POINT in_gpoint;
	POINT3D in_point3d;

	POSTGIS_DEBUGF(3, "tree1_type=%d", tree1_type);

	/* If the tree'ed argument is a polygon, do the P-i-P using the tree-based P-i-P */
	if ( tree1_type == POLYGONTYPE || tree1_type == MULTIPOLYGONTYPE )
	{
		POSTGIS_DEBUG(3, "tree is a polygon, using tree PiP");
		/* Need a gbox to calculate an outside point */
		if ( LW_FAILURE == gserialized_get_gbox_p(g1, &gbox1) )
		{
			LWGEOM* lwgeom1 = lwgeom_from_gserialized(g1);
			POSTGIS_DEBUG(3, "unable to read gbox from gserialized, calculating from scratch");
			lwgeom_calculate_gbox_geodetic(lwgeom1, &gbox1);
			lwgeom_free(lwgeom1);
		}

		/* Flip the candidate point into geographics */
		geographic_point_init(in_point->x, in_point->y, &in_gpoint);
		geog2cart(&in_gpoint, &in_point3d);

		/* If the candidate isn't in the tree box, it's not in the tree area */
		if ( ! gbox_contains_point3d(&gbox1, &in_point3d) )
		{
			POSTGIS_DEBUG(3, "in_point3d is not inside the tree gbox, CircTreePIP returning FALSE");
			return LW_FALSE;
		}
		/* The candidate point is in the box, so it *might* be inside the tree */
		else
		{
			POINT2D pt2d_outside; /* latlon */
			POINT2D pt2d_inside;
			pt2d_inside.x = in_point->x;
			pt2d_inside.y = in_point->y;
			/* Calculate a definitive outside point */
			if (gbox_pt_outside(&gbox1, &pt2d_outside) == LW_FAILURE)
				if (circ_tree_get_point_outside(tree1, &pt2d_outside) == LW_FAILURE)
					lwpgerror("%s: Unable to generate outside point!", __func__);

			POSTGIS_DEBUGF(3, "p2d_inside=POINT(%g %g) p2d_outside=POINT(%g %g)", pt2d_inside.x, pt2d_inside.y, pt2d_outside.x, pt2d_outside.y);
			/* Test the candidate point for strict containment */
			POSTGIS_DEBUG(3, "calling circ_tree_contains_point for PiP test");
			return circ_tree_contains_point(tree1, &pt2d_inside, &pt2d_outside, 0, NULL);
		}
	}
	else
	{
		POSTGIS_DEBUG(3, "tree1 not polygonal, so CircTreePIP returning FALSE");
		return LW_FALSE;
	}
}

static int
geography_distance_cache_tolerance(FunctionCallInfo fcinfo,
				   SHARED_GSERIALIZED *shared_g1,
				   SHARED_GSERIALIZED *shared_g2,
				   const SPHEROID *s,
				   double tolerance,
				   double *distance)
{
	const GSERIALIZED *g1 = shared_gserialized_get(shared_g1);
	const GSERIALIZED *g2 = shared_gserialized_get(shared_g2);
	CircTreeGeomCache* tree_cache = NULL;

	int type1 = gserialized_get_type(g1);
	int type2 = gserialized_get_type(g2);

	Assert(distance);

	/* Two points? Get outa here... */
	if ( type1 == POINTTYPE && type2 == POINTTYPE )
		return LW_FAILURE;

	/* Fetch/build our cache, if appropriate, etc... */
	tree_cache = GetCircTreeGeomCache(fcinfo, shared_g1, shared_g2);

	/* OK, we have an index at the ready! Use it for the one tree argument and */
	/* fill in the other tree argument */
	if ( tree_cache && tree_cache->gcache.argnum && tree_cache->index )
	{
		CIRC_NODE* circtree_cached = tree_cache->index;
		CIRC_NODE* circtree = NULL;
		const GSERIALIZED* g_cached;
		const GSERIALIZED* g;
		LWGEOM* lwgeom = NULL;
		int geomtype_cached;
		int geomtype;
		POINT4D p4d;

		/* We need to dynamically build a tree for the uncached side of the function call */
		if ( tree_cache->gcache.argnum == 1 )
		{
			g_cached = g1;
			g = g2;
			geomtype_cached = type1;
			geomtype = type2;
		}
		else if ( tree_cache->gcache.argnum == 2 )
		{
			g_cached = g2;
			g = g1;
			geomtype_cached = type2;
			geomtype = type1;
		}
		else
		{
			lwpgerror("geography_distance_cache this cannot happen!");
			return LW_FAILURE;
		}

		lwgeom = lwgeom_from_gserialized(g);
		if ( geomtype_cached == POLYGONTYPE || geomtype_cached == MULTIPOLYGONTYPE )
		{
			lwgeom_startpoint(lwgeom, &p4d);
			if ( CircTreePIP(circtree_cached, g_cached, &p4d) )
			{
				*distance = 0.0;
				lwgeom_free(lwgeom);
				return LW_SUCCESS;
			}
		}

		circtree = lwgeom_calculate_circ_tree(lwgeom);
		if ( geomtype == POLYGONTYPE || geomtype == MULTIPOLYGONTYPE )
		{
			POINT2D p2d;
			circ_tree_get_point(circtree_cached, &p2d);
			p4d.x = p2d.x;
			p4d.y = p2d.y;
			if ( CircTreePIP(circtree, g, &p4d) )
			{
				*distance = 0.0;
				circ_tree_free(circtree);
				lwgeom_free(lwgeom);
				return LW_SUCCESS;
			}
		}

		*distance = circ_tree_distance_tree(circtree_cached, circtree, s, tolerance);
		circ_tree_free(circtree);
		lwgeom_free(lwgeom);
		return LW_SUCCESS;
	}
	else
	{
		return LW_FAILURE;
	}
}

int
geography_distance_cache(FunctionCallInfo fcinfo,
			 SHARED_GSERIALIZED *g1,
			 SHARED_GSERIALIZED *g2,
			 const SPHEROID *s,
			 double *distance)
{
	return geography_distance_cache_tolerance(fcinfo, g1, g2, s, FP_TOLERANCE, distance);
}

int
geography_dwithin_cache(FunctionCallInfo fcinfo,
			SHARED_GSERIALIZED *g1,
			SHARED_GSERIALIZED *g2,
			const SPHEROID *s,
			double tolerance,
			int *dwithin)
{
	double distance;
	/* Ticket #2422, difference between sphere and spheroid distance can trip up the */
	/* threshold shortcircuit (stopping a calculation before the spheroid distance is actually */
	/* below the threshold. Lower in the code line, we actually reduce the threshold a little to */
	/* avoid this. */
	/* Correct fix: propogate the spheroid information all the way to the bottom of the calculation */
	/* so the "right thing" can be done in all cases. */
	if ( LW_SUCCESS == geography_distance_cache_tolerance(fcinfo, g1, g2, s, tolerance, &distance) )
	{
		*dwithin = (distance <= (tolerance + FP_TOLERANCE) ? LW_TRUE : LW_FALSE);
		return LW_SUCCESS;
	}
	return LW_FAILURE;
}

int
geography_tree_distance(const GSERIALIZED* g1, const GSERIALIZED* g2, const SPHEROID* s, double tolerance, double* distance)
{
	CIRC_NODE* circ_tree1 = NULL;
	CIRC_NODE* circ_tree2 = NULL;
	LWGEOM* lwgeom1 = NULL;
	LWGEOM* lwgeom2 = NULL;
	POINT4D pt1, pt2;

	lwgeom1 = lwgeom_from_gserialized(g1);
	lwgeom2 = lwgeom_from_gserialized(g2);
	circ_tree1 = lwgeom_calculate_circ_tree(lwgeom1);
	circ_tree2 = lwgeom_calculate_circ_tree(lwgeom2);
	lwgeom_startpoint(lwgeom1, &pt1);
	lwgeom_startpoint(lwgeom2, &pt2);

	if ( CircTreePIP(circ_tree1, g1, &pt2) || CircTreePIP(circ_tree2, g2, &pt1) )
	{
		*distance = 0.0;
	}
	else
	{
		/* Calculate tree/tree distance */
		*distance = circ_tree_distance_tree(circ_tree1, circ_tree2, s, tolerance);
	}

	circ_tree_free(circ_tree1);
	circ_tree_free(circ_tree2);
	lwgeom_free(lwgeom1);
	lwgeom_free(lwgeom2);
	return LW_SUCCESS;
}

static inline int
circ_node_is_leaf(const CIRC_NODE* node)
{
	return (node->num_nodes == 0);
}

static double
circ_node_min_distance(const CIRC_NODE* n1, const CIRC_NODE* n2)
{
	double d = sphere_distance(&(n1->center), &(n2->center));
	double r1 = n1->radius;
	double r2 = n2->radius;

	if ( d < r1 + r2 )
		return 0.0;

	return d - r1 - r2;
}

static double
circ_node_max_distance(const CIRC_NODE *n1, const CIRC_NODE *n2)
{
	return sphere_distance(&(n1->center), &(n2->center)) + n1->radius + n2->radius;
}

struct sort_node {
	CIRC_NODE *node;
	double d;
};

static int
circ_nodes_sort_cmp(const void *a, const void *b)
{
	struct sort_node *node_a = (struct sort_node *)(a);
	struct sort_node *node_b = (struct sort_node *)(b);
	if (node_a->d < node_b->d) return -1;
	else if (node_a->d > node_b->d) return 1;
	else return 0;
}

static void
circ_internal_nodes_sort(CIRC_NODE **nodes, uint32_t num_nodes, const CIRC_NODE *target_node)
{
	uint32_t i;
	struct sort_node sort_nodes[CIRC_NODE_SIZE];

	/* Copy incoming nodes into sorting array and calculate */
	/* distance to the target node */
	for (i = 0; i < num_nodes; i++)
	{
		sort_nodes[i].node = nodes[i];
		sort_nodes[i].d = sphere_distance(&(nodes[i]->center), &(target_node->center));
	}

	/* Sort the nodes and copy the result back into the input array */
	qsort(sort_nodes, num_nodes, sizeof(struct sort_node), circ_nodes_sort_cmp);
	for (i = 0; i < num_nodes; i++)
	{
		nodes[i] = sort_nodes[i].node;
	}
}


double
circ_tree_distance_tree_internal(const CIRC_NODE* n1, const CIRC_NODE* n2, double threshold,
								 double* min_dist, double* max_dist, GEOGRAPHIC_POINT* closest1, GEOGRAPHIC_POINT* closest2)
{
	double max;
	double d, d_min;
	uint32_t i;

	/* Short circuit if we've already hit the minimum */
	if( *min_dist < threshold || *min_dist == 0.0 )
		return *min_dist;

	/* If your minimum is greater than anyone's maximum, you can't hold the winner */
	if( circ_node_min_distance(n1, n2) > *max_dist )
	{
		return FLT_MAX;
	}

	/* If your maximum is a new low, we'll use that as our new global tolerance */
	max = circ_node_max_distance(n1, n2);
	if( max < *max_dist )
		*max_dist = max;

	/* Polygon on one side, primitive type on the other. Check for point-in-polygon */
	/* short circuit. */
	if ( n1->geom_type == POLYGONTYPE && n2->geom_type && ! lwtype_is_collection((uint8_t) (n2->geom_type)) )
	{
		POINT2D pt;
		circ_tree_get_point(n2, &pt);
		if ( circ_tree_contains_point(n1, &pt, &(n1->pt_outside), 0, NULL) )
		{
			*min_dist = 0.0;
			geographic_point_init(pt.x, pt.y, closest1);
			geographic_point_init(pt.x, pt.y, closest2);
			return *min_dist;
		}
	}
	/* Polygon on one side, primitive type on the other. Check for point-in-polygon */
	/* short circuit. */
	if ( n2->geom_type == POLYGONTYPE && n1->geom_type && ! lwtype_is_collection((uint8_t) (n1->geom_type)) )
	{
		POINT2D pt;
		circ_tree_get_point(n1, &pt);
		if ( circ_tree_contains_point(n2, &pt, &(n2->pt_outside), 0, NULL) )
		{
			geographic_point_init(pt.x, pt.y, closest1);
			geographic_point_init(pt.x, pt.y, closest2);
			*min_dist = 0.0;
			return *min_dist;
		}
	}

	/* Both leaf nodes, do a real distance calculation */
	if( circ_node_is_leaf(n1) && circ_node_is_leaf(n2) )
	{
		double d;
		GEOGRAPHIC_POINT close1, close2;
		/* One of the nodes is a point */
		if ( n1->p1 == n1->p2 || n2->p1 == n2->p2 )
		{
			GEOGRAPHIC_EDGE e;
			GEOGRAPHIC_POINT gp1, gp2;

			/* Both nodes are points! */
			if ( n1->p1 == n1->p2 && n2->p1 == n2->p2 )
			{
				geographic_point_init(n1->p1->x, n1->p1->y, &gp1);
				geographic_point_init(n2->p1->x, n2->p1->y, &gp2);
				close1 = gp1; close2 = gp2;
				d = sphere_distance(&gp1, &gp2);
			}
				/* Node 1 is a point */
			else if ( n1->p1 == n1->p2 )
			{
				geographic_point_init(n1->p1->x, n1->p1->y, &gp1);
				geographic_point_init(n2->p1->x, n2->p1->y, &(e.start));
				geographic_point_init(n2->p2->x, n2->p2->y, &(e.end));
				close1 = gp1;
				d = edge_distance_to_point(&e, &gp1, &close2);
			}
				/* Node 2 is a point */
			else
			{
				/* FIX
                geographic_point_init(n2->p1->x, n2->p1->y, &gp1); */
				geographic_point_init(n2->p1->x, n2->p1->y, &gp2);
				geographic_point_init(n1->p1->x, n1->p1->y, &(e.start));
				geographic_point_init(n1->p2->x, n1->p2->y, &(e.end));
				/* FIX
                close1 = gp1;
                d = edge_distance_to_point(&e, &gp1, &close2); */
				close2 = gp2;
				d = edge_distance_to_point(&e, &gp2, &close1);
			}
		}
			/* Both nodes are edges */
		else
		{
			GEOGRAPHIC_EDGE e1, e2;
			GEOGRAPHIC_POINT g;
			POINT3D A1, A2, B1, B2;
			geographic_point_init(n1->p1->x, n1->p1->y, &(e1.start));
			geographic_point_init(n1->p2->x, n1->p2->y, &(e1.end));
			geographic_point_init(n2->p1->x, n2->p1->y, &(e2.start));
			geographic_point_init(n2->p2->x, n2->p2->y, &(e2.end));
			geog2cart(&(e1.start), &A1);
			geog2cart(&(e1.end), &A2);
			geog2cart(&(e2.start), &B1);
			geog2cart(&(e2.end), &B2);
			if ( edge_intersects(&A1, &A2, &B1, &B2) )
			{
				d = 0.0;
				edge_intersection(&e1, &e2, &g);
				close1 = close2 = g;
			}
			else
			{
				d = edge_distance_to_edge(&e1, &e2, &close1, &close2);
			}
		}
		if ( d < *min_dist )
		{
			*min_dist = d;
			*closest1 = close1;
			*closest2 = close2;
		}
		return d;
	}
	else
	{
		d_min = FLT_MAX;
		/* Drive the recursion into the COLLECTION types first so we end up with */
		/* pairings of primitive geometries that can be forced into the point-in-polygon */
		/* tests above. */
		if ( n1->geom_type && lwtype_is_collection((uint8_t) (n1->geom_type)) )
		{
			circ_internal_nodes_sort(n1->nodes, n1->num_nodes, n2);
			for ( i = 0; i < n1->num_nodes; i++ )
			{
				d = circ_tree_distance_tree_internal(n1->nodes[i], n2, threshold, min_dist, max_dist, closest1, closest2);
				d_min = FP_MIN(d_min, d);
			}
		}
		else if ( n2->geom_type && lwtype_is_collection((uint8_t) (n2->geom_type)) )
		{
			circ_internal_nodes_sort(n2->nodes, n2->num_nodes, n1);
			for ( i = 0; i < n2->num_nodes; i++ )
			{
				d = circ_tree_distance_tree_internal(n1, n2->nodes[i], threshold, min_dist, max_dist, closest1, closest2);
				d_min = FP_MIN(d_min, d);
			}
		}
		else if ( ! circ_node_is_leaf(n1) )
		{
			circ_internal_nodes_sort(n1->nodes, n1->num_nodes, n2);
			for ( i = 0; i < n1->num_nodes; i++ )
			{
				d = circ_tree_distance_tree_internal(n1->nodes[i], n2, threshold, min_dist, max_dist, closest1, closest2);
				d_min = FP_MIN(d_min, d);
			}
		}
		else if ( ! circ_node_is_leaf(n2) )
		{
			circ_internal_nodes_sort(n2->nodes, n2->num_nodes, n1);
			for ( i = 0; i < n2->num_nodes; i++ )
			{
				d = circ_tree_distance_tree_internal(n1, n2->nodes[i], threshold, min_dist, max_dist, closest1, closest2);
				d_min = FP_MIN(d_min, d);
			}
		}
		else
		{
			/* Never get here */
		}

		return d_min;
	}
}

/* Closest point for geographies. */
LWGEOM *
geography_tree_closestpoint(const GSERIALIZED* g1, const GSERIALIZED* g2, double threshold)
{
    CIRC_NODE* circ_tree1 = NULL;
    CIRC_NODE* circ_tree2 = NULL;
    LWGEOM* lwgeom1 = NULL;
    LWGEOM* lwgeom2 = NULL;
    double min_dist = FLT_MAX;
    double max_dist = FLT_MAX;
    GEOGRAPHIC_POINT closest1, closest2;
    LWGEOM *result;
    POINT4D p;

    lwgeom1 = lwgeom_from_gserialized(g1);
    lwgeom2 = lwgeom_from_gserialized(g2);
    circ_tree1 = lwgeom_calculate_circ_tree(lwgeom1);
    circ_tree2 = lwgeom_calculate_circ_tree(lwgeom2);

    circ_tree_distance_tree_internal(circ_tree1, circ_tree2, threshold,
                                     &min_dist, &max_dist, &closest1, &closest2);

    p.x = rad2deg(closest1.lon);
    p.y = rad2deg(closest1.lat);
    result = (LWGEOM *)lwpoint_make2d(gserialized_get_srid(g1), p.x, p.y);

    circ_tree_free(circ_tree1);
    circ_tree_free(circ_tree2);
    lwgeom_free(lwgeom1);
    lwgeom_free(lwgeom2);
    return result;
}

/* Shortest line for geographies. */
LWGEOM *
geography_tree_shortestline(const GSERIALIZED* g1, const GSERIALIZED* g2, double threshold)
{
    CIRC_NODE* circ_tree1 = NULL;
    CIRC_NODE* circ_tree2 = NULL;
    LWGEOM* lwgeom1 = NULL;
    LWGEOM* lwgeom2 = NULL;
    double min_dist = FLT_MAX;
    double max_dist = FLT_MAX;
    GEOGRAPHIC_POINT closest1, closest2;
    LWGEOM *geoms[2];
    LWGEOM *result;
    POINT4D p1, p2;

    lwgeom1 = lwgeom_from_gserialized(g1);
    lwgeom2 = lwgeom_from_gserialized(g2);
    circ_tree1 = lwgeom_calculate_circ_tree(lwgeom1);
    circ_tree2 = lwgeom_calculate_circ_tree(lwgeom2);

    circ_tree_distance_tree_internal(circ_tree1, circ_tree2, threshold,
                                     &min_dist, &max_dist, &closest1, &closest2);

    p1.x = rad2deg(closest1.lon);
    p1.y = rad2deg(closest1.lat);
    p2.x = rad2deg(closest2.lon);
    p2.y = rad2deg(closest2.lat);

    geoms[0] = (LWGEOM *)lwpoint_make2d(gserialized_get_srid(g1), p1.x, p1.y);
    geoms[1] = (LWGEOM *)lwpoint_make2d(gserialized_get_srid(g1), p2.x, p2.y);
    result = (LWGEOM *)lwline_from_lwgeom_array(geoms[0]->srid, 2, geoms);

    lwgeom_free(geoms[0]);
    lwgeom_free(geoms[1]);
    circ_tree_free(circ_tree1);
    circ_tree_free(circ_tree2);
    lwgeom_free(lwgeom1);
    lwgeom_free(lwgeom2);
    return result;
}
