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
 * Copyright (C) 2009 Paul Ramsey <pramsey@cleverelephant.ca>
 *
 **********************************************************************/


#include "postgres.h"
#include "catalog/pg_type.h" /* for CSTRINGOID */

#include "../postgis_config.h"

#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>

#include "liblwgeom.h"                  /* For standard geometry types. */
#include "liblwgeom_internal.h"         /* For FP comparators. */
#include "lwgeom_pg.h"                  /* For debugging macros. */
#include "geography.h"                  /* For utility functions. */
#include "geography_measurement_trees.h" /* For circ_tree caching */
#include "lwgeom_transform.h"            /* For SRID functions */

#ifdef PROJ_GEODESIC
/* round to 10 nm precision */
#define INVMINDIST 1.0e8
#else
/* round to 100 nm precision */
#define INVMINDIST 1.0e9
#endif

Datum geography_distance(PG_FUNCTION_ARGS);
Datum geography_distance_uncached(PG_FUNCTION_ARGS);
Datum geography_distance_knn(PG_FUNCTION_ARGS);
Datum geography_distance_tree(PG_FUNCTION_ARGS);
Datum geography_dwithin(PG_FUNCTION_ARGS);
Datum geography_dwithin_uncached(PG_FUNCTION_ARGS);
Datum geography_area(PG_FUNCTION_ARGS);
Datum geography_length(PG_FUNCTION_ARGS);
Datum geography_expand(PG_FUNCTION_ARGS);
Datum geography_point_outside(PG_FUNCTION_ARGS);
Datum geography_covers(PG_FUNCTION_ARGS);
Datum geography_coveredby(PG_FUNCTION_ARGS);
Datum geography_bestsrid(PG_FUNCTION_ARGS);
Datum geography_perimeter(PG_FUNCTION_ARGS);
Datum geography_project(PG_FUNCTION_ARGS);
Datum geography_azimuth(PG_FUNCTION_ARGS);
Datum geography_segmentize(PG_FUNCTION_ARGS);
Datum geography_line_interpolate_point(PG_FUNCTION_ARGS);
Datum geography_line_locate_point(PG_FUNCTION_ARGS);
Datum geography_closestpoint(PG_FUNCTION_ARGS);
Datum geography_shortestline(PG_FUNCTION_ARGS);
Datum geography_segmentize1(PG_FUNCTION_ARGS);


PG_FUNCTION_INFO_V1(geography_distance_knn);
Datum geography_distance_knn(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom1 = NULL;
	LWGEOM *lwgeom2 = NULL;
	GSERIALIZED *g1 = NULL;
	GSERIALIZED *g2 = NULL;
	double distance;
	double tolerance = FP_TOLERANCE;
	bool use_spheroid = false; /* must use sphere, can't get index to harmonize with spheroid */
	SPHEROID s;

	/* Get our geometry objects loaded into memory. */
	g1 = PG_GETARG_GSERIALIZED_P(0);
	g2 = PG_GETARG_GSERIALIZED_P(1);

	gserialized_error_if_srid_mismatch(g1, g2, __func__);

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g1), &s);

	/* Set to sphere if requested */
	if ( ! use_spheroid )
		s.a = s.b = s.radius;

	lwgeom1 = lwgeom_from_gserialized(g1);
	lwgeom2 = lwgeom_from_gserialized(g2);

	/* Return NULL on empty arguments. */
	if ( lwgeom_is_empty(lwgeom1) || lwgeom_is_empty(lwgeom2) )
	{
		PG_FREE_IF_COPY(g1, 0);
		PG_FREE_IF_COPY(g2, 1);
		PG_RETURN_NULL();
	}

	/* Make sure we have boxes attached */
	lwgeom_add_bbox_deep(lwgeom1, NULL);
	lwgeom_add_bbox_deep(lwgeom2, NULL);

	distance = lwgeom_distance_spheroid(lwgeom1, lwgeom2, &s, tolerance);

	POSTGIS_DEBUGF(2, "[GIST] '%s' got distance %g", __func__, distance);

	/* Clean up */
	lwgeom_free(lwgeom1);
	lwgeom_free(lwgeom2);
	PG_FREE_IF_COPY(g1, 0);
	PG_FREE_IF_COPY(g2, 1);

	/* Something went wrong, negative return... should already be eloged, return NULL */
	if ( distance < 0.0 )
	{
		PG_RETURN_NULL();
	}

	PG_RETURN_FLOAT8(distance);
}

/*
** geography_distance_uncached(GSERIALIZED *g1, GSERIALIZED *g2, double tolerance, boolean use_spheroid)
** returns double distance in meters
*/
PG_FUNCTION_INFO_V1(geography_distance_uncached);
Datum geography_distance_uncached(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom1 = NULL;
	LWGEOM *lwgeom2 = NULL;
	GSERIALIZED *g1 = NULL;
	GSERIALIZED *g2 = NULL;
	double distance;
	double tolerance = FP_TOLERANCE;
	bool use_spheroid = true;
	SPHEROID s;

	/* Get our geometry objects loaded into memory. */
	g1 = PG_GETARG_GSERIALIZED_P(0);
	g2 = PG_GETARG_GSERIALIZED_P(1);

	/* Read our tolerance value. */
	if ( PG_NARGS() > 2 && ! PG_ARGISNULL(2) )
		tolerance = PG_GETARG_FLOAT8(2);

	/* Read our calculation type. */
	if ( PG_NARGS() > 3 && ! PG_ARGISNULL(3) )
		use_spheroid = PG_GETARG_BOOL(3);

	gserialized_error_if_srid_mismatch(g1, g2, __func__);

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g1), &s);

	/* Set to sphere if requested */
	if ( ! use_spheroid )
		s.a = s.b = s.radius;

	lwgeom1 = lwgeom_from_gserialized(g1);
	lwgeom2 = lwgeom_from_gserialized(g2);

	/* Return NULL on empty arguments. */
	if ( lwgeom_is_empty(lwgeom1) || lwgeom_is_empty(lwgeom2) )
	{
		PG_FREE_IF_COPY(g1, 0);
		PG_FREE_IF_COPY(g2, 1);
		PG_RETURN_NULL();
	}

	/* Make sure we have boxes attached */
	lwgeom_add_bbox_deep(lwgeom1, NULL);
	lwgeom_add_bbox_deep(lwgeom2, NULL);

	distance = lwgeom_distance_spheroid(lwgeom1, lwgeom2, &s, tolerance);

	POSTGIS_DEBUGF(2, "[GIST] '%s' got distance %g", __func__, distance);

	/* Clean up */
	lwgeom_free(lwgeom1);
	lwgeom_free(lwgeom2);
	PG_FREE_IF_COPY(g1, 0);
	PG_FREE_IF_COPY(g2, 1);

	/* Something went wrong, negative return... should already be eloged, return NULL */
	if ( distance < 0.0 )
	{
		PG_RETURN_NULL();
	}

	PG_RETURN_FLOAT8(distance);
}


/*
** geography_distance(GSERIALIZED *g1, GSERIALIZED *g2, double tolerance, boolean use_spheroid)
** returns double distance in meters
*/
PG_FUNCTION_INFO_V1(geography_distance);
Datum geography_distance(PG_FUNCTION_ARGS)
{
	SHARED_GSERIALIZED *shared_geom1 = ToastCacheGetGeometry(fcinfo, 0);
	SHARED_GSERIALIZED *shared_geom2 = ToastCacheGetGeometry(fcinfo, 1);
	const GSERIALIZED *g1 = shared_gserialized_get(shared_geom1);
	const GSERIALIZED *g2 = shared_gserialized_get(shared_geom2);
	double distance;
	bool use_spheroid = true;
	SPHEROID s;


	if (PG_NARGS() > 2)
		use_spheroid = PG_GETARG_BOOL(2);

	gserialized_error_if_srid_mismatch(g1, g2, __func__);

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g1), &s);

	/* Set to sphere if requested */
	if ( ! use_spheroid )
		s.a = s.b = s.radius;

	/* Return NULL on empty arguments. */
	if ( gserialized_is_empty(g1) || gserialized_is_empty(g2) )
	{
		PG_RETURN_NULL();
	}

	/* Do the brute force calculation if the cached calculation doesn't tick over */
	if (LW_FAILURE == geography_distance_cache(fcinfo, shared_geom1, shared_geom2, &s, &distance))
	{
		/* default to using tree-based distance calculation at all times */
		/* in standard distance call. */
		geography_tree_distance(g1, g2, &s, FP_TOLERANCE, &distance);
		/*
		LWGEOM* lwgeom1 = lwgeom_from_gserialized(g1);
		LWGEOM* lwgeom2 = lwgeom_from_gserialized(g2);
		distance = lwgeom_distance_spheroid(lwgeom1, lwgeom2, &s, tolerance);
		lwgeom_free(lwgeom1);
		lwgeom_free(lwgeom2);
		*/
	}

	/* Knock off any funny business at the nanometer level, ticket #2168 */
	distance = round(distance * INVMINDIST) / INVMINDIST;

	/* Something went wrong, negative return... should already be eloged, return NULL */
	if ( distance < 0.0 )
	{
		elog(ERROR, "distance returned negative!");
		PG_RETURN_NULL();
	}

	PG_RETURN_FLOAT8(distance);
}

/*
** geography_dwithin(GSERIALIZED *g1, GSERIALIZED *g2, double tolerance, boolean use_spheroid)
** returns double distance in meters
*/
PG_FUNCTION_INFO_V1(geography_dwithin);
Datum geography_dwithin(PG_FUNCTION_ARGS)
{
	SHARED_GSERIALIZED *shared_geom1 = ToastCacheGetGeometry(fcinfo, 0);
	SHARED_GSERIALIZED *shared_geom2 = ToastCacheGetGeometry(fcinfo, 1);
	const GSERIALIZED *g1 = shared_gserialized_get(shared_geom1);
	const GSERIALIZED *g2 = shared_gserialized_get(shared_geom2);
	SPHEROID s;
	double tolerance = 0.0;
	bool use_spheroid = true;
	double distance;
	int dwithin = LW_FALSE;

	gserialized_error_if_srid_mismatch(g1, g2, __func__);

	/* Read our tolerance value. */
	if ( PG_NARGS() > 2 && ! PG_ARGISNULL(2) )
		tolerance = PG_GETARG_FLOAT8(2);

	/* Read our calculation type. */
	if ( PG_NARGS() > 3 && ! PG_ARGISNULL(3) )
		use_spheroid = PG_GETARG_BOOL(3);

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g1), &s);

	/* Set to sphere if requested */
	if ( ! use_spheroid )
		s.a = s.b = s.radius;

	/* Return FALSE on empty arguments. */
	if ( gserialized_is_empty(g1) || gserialized_is_empty(g2) )
		PG_RETURN_BOOL(false);

	/* Do the brute force calculation if the cached calculation doesn't tick over */
	if (LW_FAILURE == geography_dwithin_cache(fcinfo, shared_geom1, shared_geom2, &s, tolerance, &dwithin))
	{
		LWGEOM* lwgeom1 = lwgeom_from_gserialized(g1);
		LWGEOM* lwgeom2 = lwgeom_from_gserialized(g2);
		distance = lwgeom_distance_spheroid(lwgeom1, lwgeom2, &s, tolerance);
		/* Something went wrong... */
		if ( distance < 0.0 )
			elog(ERROR, "lwgeom_distance_spheroid returned negative!");
		dwithin = (distance <= tolerance);
		lwgeom_free(lwgeom1);
		lwgeom_free(lwgeom2);
	}

	PG_RETURN_BOOL(dwithin);
}

PG_FUNCTION_INFO_V1(geography_intersects);
Datum geography_intersects(PG_FUNCTION_ARGS)
{
	PG_RETURN_BOOL(CallerFInfoFunctionCall2(
		geography_dwithin, fcinfo->flinfo, InvalidOid,
		PG_GETARG_DATUM(0), PG_GETARG_DATUM(1)));
}

/*
** geography_distance_tree(GSERIALIZED *g1, GSERIALIZED *g2, double tolerance, boolean use_spheroid)
** returns double distance in meters
*/
PG_FUNCTION_INFO_V1(geography_distance_tree);
Datum geography_distance_tree(PG_FUNCTION_ARGS)
{
	GSERIALIZED *g1 = NULL;
	GSERIALIZED *g2 = NULL;
	double tolerance = 0.0;
	double distance;
	bool use_spheroid = true;
	SPHEROID s;

	/* Get our geometry objects loaded into memory. */
	g1 = PG_GETARG_GSERIALIZED_P(0);
	g2 = PG_GETARG_GSERIALIZED_P(1);

	gserialized_error_if_srid_mismatch(g1, g2, __func__);

	/* Return zero on empty arguments. */
	if ( gserialized_is_empty(g1) || gserialized_is_empty(g2) )
	{
		PG_FREE_IF_COPY(g1, 0);
		PG_FREE_IF_COPY(g2, 1);
		PG_RETURN_FLOAT8(0.0);
	}

	/* Read our tolerance value. */
	if ( PG_NARGS() > 2 && ! PG_ARGISNULL(2) )
		tolerance = PG_GETARG_FLOAT8(2);

	/* Read our calculation type. */
	if ( PG_NARGS() > 3 && ! PG_ARGISNULL(3) )
		use_spheroid = PG_GETARG_BOOL(3);

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g1), &s);

	/* Set to sphere if requested */
	if ( ! use_spheroid )
		s.a = s.b = s.radius;

	if  ( geography_tree_distance(g1, g2, &s, tolerance, &distance) == LW_FAILURE )
	{
		elog(ERROR, "geography_distance_tree failed!");
		PG_RETURN_NULL();
	}

	PG_RETURN_FLOAT8(distance);
}



/*
** geography_dwithin_uncached(GSERIALIZED *g1, GSERIALIZED *g2, double tolerance, boolean use_spheroid)
** returns double distance in meters
*/
PG_FUNCTION_INFO_V1(geography_dwithin_uncached);
Datum geography_dwithin_uncached(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom1 = NULL;
	LWGEOM *lwgeom2 = NULL;
	GSERIALIZED *g1 = NULL;
	GSERIALIZED *g2 = NULL;
	double tolerance = 0.0;
	double distance;
	bool use_spheroid = true;
	SPHEROID s;

	/* Get our geometry objects loaded into memory. */
	g1 = PG_GETARG_GSERIALIZED_P(0);
	g2 = PG_GETARG_GSERIALIZED_P(1);
	gserialized_error_if_srid_mismatch(g1, g2, __func__);

	/* Read our tolerance value. */
	if ( PG_NARGS() > 2 && ! PG_ARGISNULL(2) )
		tolerance = PG_GETARG_FLOAT8(2);

	/* Read our calculation type. */
	if ( PG_NARGS() > 3 && ! PG_ARGISNULL(3) )
		use_spheroid = PG_GETARG_BOOL(3);

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g1), &s);

	/* Set to sphere if requested */
	if ( ! use_spheroid )
		s.a = s.b = s.radius;

	lwgeom1 = lwgeom_from_gserialized(g1);
	lwgeom2 = lwgeom_from_gserialized(g2);

	/* Return FALSE on empty arguments. */
	if ( lwgeom_is_empty(lwgeom1) || lwgeom_is_empty(lwgeom2) )
	{
		PG_RETURN_BOOL(false);
	}

	distance = lwgeom_distance_spheroid(lwgeom1, lwgeom2, &s, tolerance);

	/* Clean up */
	lwgeom_free(lwgeom1);
	lwgeom_free(lwgeom2);
	PG_FREE_IF_COPY(g1, 0);
	PG_FREE_IF_COPY(g2, 1);

	/* Something went wrong... should already be eloged, return FALSE */
	if ( distance < 0.0 )
	{
		elog(ERROR, "lwgeom_distance_spheroid returned negative!");
		PG_RETURN_BOOL(false);
	}

	PG_RETURN_BOOL(distance <= tolerance);
}


/*
** geography_expand(GSERIALIZED *g) returns *GSERIALIZED
**
** warning, this tricky little function does not expand the
** geometry at all, just re-writes bounding box value to be
** a bit bigger. only useful when passing the result along to
** an index operator (&&)
*/
PG_FUNCTION_INFO_V1(geography_expand);
Datum geography_expand(PG_FUNCTION_ARGS)
{
	GSERIALIZED *g = NULL;
	GSERIALIZED *g_out = NULL;
	double unit_distance, distance;

	/* Get a wholly-owned pointer to the geography */
	g = PG_GETARG_GSERIALIZED_P_COPY(0);

	/* Read our distance value and normalize to unit-sphere. */
	distance = PG_GETARG_FLOAT8(1);
	/* Magic 1% expansion is to bridge difference between potential */
	/* spheroidal input distance and fact that expanded box filter is */
	/* calculated on sphere */
	unit_distance = 1.01 * distance / WGS84_RADIUS;

	/* Try the expansion */
	g_out = gserialized_expand(g, unit_distance);

	/* If the expansion fails, the return our input */
	if ( g_out == NULL )
	{
		PG_RETURN_POINTER(g);
	}

	if ( g_out != g )
	{
		pfree(g);
	}

	PG_RETURN_POINTER(g_out);
}

/*
** geography_area(GSERIALIZED *g)
** returns double area in meters square
*/
PG_FUNCTION_INFO_V1(geography_area);
Datum geography_area(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom = NULL;
	GSERIALIZED *g = NULL;
	GBOX gbox;
	double area;
	bool use_spheroid = LW_TRUE;
	SPHEROID s;

	/* Get our geometry object loaded into memory. */
	g = PG_GETARG_GSERIALIZED_P(0);

	/* Read our calculation type */
	use_spheroid = PG_GETARG_BOOL(1);

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g), &s);

	lwgeom = lwgeom_from_gserialized(g);

	/* EMPTY things have no area */
	if ( lwgeom_is_empty(lwgeom) )
	{
		lwgeom_free(lwgeom);
		PG_RETURN_FLOAT8(0.0);
	}

	if ( lwgeom->bbox )
		gbox = *(lwgeom->bbox);
	else
		lwgeom_calculate_gbox_geodetic(lwgeom, &gbox);

#ifndef PROJ_GEODESIC
	/* Test for cases that are currently not handled by spheroid code */
	if ( use_spheroid )
	{
		/* We can't circle the poles right now */
		if ( FP_GTEQ(gbox.zmax,1.0) || FP_LTEQ(gbox.zmin,-1.0) )
			use_spheroid = LW_FALSE;
		/* We can't cross the equator right now */
		if ( gbox.zmax > 0.0 && gbox.zmin < 0.0 )
			use_spheroid = LW_FALSE;
	}
#endif /* ifndef PROJ_GEODESIC */

	/* User requests spherical calculation, turn our spheroid into a sphere */
	if ( ! use_spheroid )
		s.a = s.b = s.radius;

	/* Calculate the area */
	if ( use_spheroid )
		area = lwgeom_area_spheroid(lwgeom, &s);
	else
		area = lwgeom_area_sphere(lwgeom, &s);

	/* Clean up */
	lwgeom_free(lwgeom);
	PG_FREE_IF_COPY(g, 0);

	/* Something went wrong... */
	if ( area < 0.0 )
	{
		elog(ERROR, "lwgeom_area_spher(oid) returned area < 0.0");
		PG_RETURN_NULL();
	}

	PG_RETURN_FLOAT8(area);
}

/*
** geography_perimeter(GSERIALIZED *g)
** returns double perimeter in meters for area features
*/
PG_FUNCTION_INFO_V1(geography_perimeter);
Datum geography_perimeter(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom = NULL;
	GSERIALIZED *g = NULL;
	double length;
	bool use_spheroid = LW_TRUE;
	SPHEROID s;
	int type;

	/* Get our geometry object loaded into memory. */
	g = PG_GETARG_GSERIALIZED_P(0);

	/* Only return for area features. */
	type = gserialized_get_type(g);
	if ( ! (type == POLYGONTYPE || type == MULTIPOLYGONTYPE || type == COLLECTIONTYPE) )
	{
		PG_RETURN_FLOAT8(0.0);
	}

	lwgeom = lwgeom_from_gserialized(g);

	/* EMPTY things have no perimeter */
	if ( lwgeom_is_empty(lwgeom) )
	{
		lwgeom_free(lwgeom);
		PG_RETURN_FLOAT8(0.0);
	}

	/* Read our calculation type */
	use_spheroid = PG_GETARG_BOOL(1);

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g), &s);

	/* User requests spherical calculation, turn our spheroid into a sphere */
	if ( ! use_spheroid )
		s.a = s.b = s.radius;

	/* Calculate the length */
	length = lwgeom_length_spheroid(lwgeom, &s);

	/* Something went wrong... */
	if ( length < 0.0 )
	{
		elog(ERROR, "lwgeom_length_spheroid returned length < 0.0");
		PG_RETURN_NULL();
	}

	/* Clean up, but not all the way to the point arrays */
	lwgeom_free(lwgeom);

	PG_FREE_IF_COPY(g, 0);
	PG_RETURN_FLOAT8(length);
}

/*
** geography_length(GSERIALIZED *g)
** returns double length in meters
*/
PG_FUNCTION_INFO_V1(geography_length);
Datum geography_length(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom = NULL;
	GSERIALIZED *g = NULL;
	double length;
	bool use_spheroid = LW_TRUE;
	SPHEROID s;

	/* Get our geometry object loaded into memory. */
	g = PG_GETARG_GSERIALIZED_P(0);
	lwgeom = lwgeom_from_gserialized(g);

	/* EMPTY things have no length */
	if ( lwgeom_is_empty(lwgeom) || lwgeom->type == POLYGONTYPE || lwgeom->type == MULTIPOLYGONTYPE )
	{
		lwgeom_free(lwgeom);
		PG_RETURN_FLOAT8(0.0);
	}

	/* Read our calculation type */
	use_spheroid = PG_GETARG_BOOL(1);

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g), &s);

	/* User requests spherical calculation, turn our spheroid into a sphere */
	if ( ! use_spheroid )
		s.a = s.b = s.radius;

	/* Calculate the length */
	length = lwgeom_length_spheroid(lwgeom, &s);

	/* Something went wrong... */
	if ( length < 0.0 )
	{
		elog(ERROR, "lwgeom_length_spheroid returned length < 0.0");
		PG_RETURN_NULL();
	}

	/* Clean up */
	lwgeom_free(lwgeom);

	PG_FREE_IF_COPY(g, 0);
	PG_RETURN_FLOAT8(length);
}

/*
** geography_point_outside(GSERIALIZED *g)
** returns point outside the object
*/
PG_FUNCTION_INFO_V1(geography_point_outside);
Datum geography_point_outside(PG_FUNCTION_ARGS)
{
	GBOX gbox;
	GSERIALIZED *g = NULL;
	GSERIALIZED *g_out = NULL;
	LWGEOM *lwpoint = NULL;
	POINT2D pt;

	/* Get our geometry object loaded into memory. */
	g = PG_GETARG_GSERIALIZED_P(0);

	/* We need the bounding box to get an outside point for area algorithm */
	if ( gserialized_get_gbox_p(g, &gbox) == LW_FAILURE )
	{
		POSTGIS_DEBUG(4,"gserialized_get_gbox_p returned LW_FAILURE");
		elog(ERROR, "Error in gserialized_get_gbox_p calculation.");
		PG_RETURN_NULL();
	}
	POSTGIS_DEBUGF(4, "got gbox %s", gbox_to_string(&gbox));

	/* Get an exterior point, based on this gbox */
	gbox_pt_outside(&gbox, &pt);

	lwpoint = (LWGEOM*) lwpoint_make2d(4326, pt.x, pt.y);

	g_out = geography_serialize(lwpoint);

	PG_FREE_IF_COPY(g, 0);
	PG_RETURN_POINTER(g_out);

}

/*
** geography_covers(GSERIALIZED *g, GSERIALIZED *g) returns boolean
** Only works for (multi)points and (multi)polygons currently.
** Attempts a simple point-in-polygon test on the polygon and point.
** Current algorithm does not distinguish between points on edge
** and points within.
*/
PG_FUNCTION_INFO_V1(geography_covers);
Datum geography_covers(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom1 = NULL;
	LWGEOM *lwgeom2 = NULL;
	GSERIALIZED *g1 = NULL;
	GSERIALIZED *g2 = NULL;
	int result = LW_FALSE;

	/* Get our geometry objects loaded into memory. */
	g1 = PG_GETARG_GSERIALIZED_P(0);
	g2 = PG_GETARG_GSERIALIZED_P(1);
	gserialized_error_if_srid_mismatch(g1, g2, __func__);

	/* Construct our working geometries */
	lwgeom1 = lwgeom_from_gserialized(g1);
	lwgeom2 = lwgeom_from_gserialized(g2);

	/* EMPTY never intersects with another geometry */
	if ( lwgeom_is_empty(lwgeom1) || lwgeom_is_empty(lwgeom2) )
	{
		lwgeom_free(lwgeom1);
		lwgeom_free(lwgeom2);
		PG_FREE_IF_COPY(g1, 0);
		PG_FREE_IF_COPY(g2, 1);
		PG_RETURN_BOOL(false);
	}

	/* Calculate answer */
	result = lwgeom_covers_lwgeom_sphere(lwgeom1, lwgeom2);

	/* Clean up */
	lwgeom_free(lwgeom1);
	lwgeom_free(lwgeom2);
	PG_FREE_IF_COPY(g1, 0);
	PG_FREE_IF_COPY(g2, 1);

	PG_RETURN_BOOL(result);
}

PG_FUNCTION_INFO_V1(geography_coveredby);
Datum geography_coveredby(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom1 = NULL;
	LWGEOM *lwgeom2 = NULL;
	GSERIALIZED *g1 = NULL;
	GSERIALIZED *g2 = NULL;
	int result = LW_FALSE;

	/* Get our geometry objects loaded into memory. */
	/* Pick them up in reverse order to covers */
	g1 = PG_GETARG_GSERIALIZED_P(1);
	g2 = PG_GETARG_GSERIALIZED_P(0);
	gserialized_error_if_srid_mismatch(g1, g2, __func__);

	/* Construct our working geometries */
	lwgeom1 = lwgeom_from_gserialized(g1);
	lwgeom2 = lwgeom_from_gserialized(g2);

	/* EMPTY never intersects with another geometry */
	if ( lwgeom_is_empty(lwgeom1) || lwgeom_is_empty(lwgeom2) )
	{
		lwgeom_free(lwgeom1);
		lwgeom_free(lwgeom2);
		PG_FREE_IF_COPY(g1, 1);
		PG_FREE_IF_COPY(g2, 0);
		PG_RETURN_BOOL(false);
	}

	/* Calculate answer */
	result = lwgeom_covers_lwgeom_sphere(lwgeom1, lwgeom2);

	/* Clean up */
	lwgeom_free(lwgeom1);
	lwgeom_free(lwgeom2);
	PG_FREE_IF_COPY(g1, 1);
	PG_FREE_IF_COPY(g2, 0);

	PG_RETURN_BOOL(result);
}

/*
** geography_bestsrid(GSERIALIZED *g, GSERIALIZED *g) returns int
** Utility function. Returns negative SRID numbers that match to the
** numbers handled in code by the transform(lwgeom, srid) function.
** UTM, polar stereographic and mercator as fallback. To be used
** in wrapping existing geometry functions in SQL to provide access
** to them in the geography module.
*/
PG_FUNCTION_INFO_V1(geography_bestsrid);
Datum geography_bestsrid(PG_FUNCTION_ARGS)
{
	GBOX gbox, gbox1, gbox2;
	GSERIALIZED *g1 = NULL;
	GSERIALIZED *g2 = NULL;
	int empty1 = LW_FALSE;
	int empty2 = LW_FALSE;
	double xwidth, ywidth;
	POINT2D center;
	Datum d1 = PG_GETARG_DATUM(0);

	/* Get our geometry objects loaded into memory. */
	g1 = (GSERIALIZED*)PG_DETOAST_DATUM(d1);
	/* Synchronize our box types */
	gbox1.flags = gserialized_get_lwflags(g1);
	/* Calculate if the geometry is empty. */
	empty1 = gserialized_is_empty(g1);
	/* Calculate a geocentric bounds for the objects */
	if ( ! empty1 && gserialized_get_gbox_p(g1, &gbox1) == LW_FAILURE )
		elog(ERROR, "Error in geography_bestsrid calling gserialized_get_gbox_p(g1, &gbox1)");

	POSTGIS_DEBUGF(4, "calculated gbox = %s", gbox_to_string(&gbox1));

	/* If we have a unique second argument, fill in all the necessary variables. */
	if (PG_NARGS() > 1)
	{
		Datum d2 = PG_GETARG_DATUM(1);
		g2 = (GSERIALIZED*)PG_DETOAST_DATUM(d2);
		gbox2.flags = gserialized_get_lwflags(g2);
		empty2 = gserialized_is_empty(g2);
		if ( ! empty2 && gserialized_get_gbox_p(g2, &gbox2) == LW_FAILURE )
			elog(ERROR, "Error in geography_bestsrid calling gserialized_get_gbox_p(g2, &gbox2)");
	}
	/*
	** If no unique second argument, copying the box for the first
	** argument will give us the right answer for all subsequent tests.
	*/
	else
	{
		gbox = gbox2 = gbox1;
	}

	/* Both empty? We don't have an answer. */
	if ( empty1 && empty2 )
		PG_RETURN_NULL();

	/* One empty? We can use the other argument values as infill. Otherwise merge the boxen */
	if ( empty1 )
		gbox = gbox2;
	else if ( empty2 )
		gbox = gbox1;
	else
		gbox_union(&gbox1, &gbox2, &gbox);

	gbox_centroid(&gbox, &center);

	/* Width and height in degrees */
	xwidth = 180.0 * gbox_angular_width(&gbox)  / M_PI;
	ywidth = 180.0 * gbox_angular_height(&gbox) / M_PI;

	POSTGIS_DEBUGF(2, "xwidth %g", xwidth);
	POSTGIS_DEBUGF(2, "ywidth %g", ywidth);
	POSTGIS_DEBUGF(2, "center POINT(%g %g)", center.x, center.y);

	/* Are these data arctic? Lambert Azimuthal Equal Area North. */
	if ( center.y > 70.0 && ywidth < 45.0 )
	{
		PG_RETURN_INT32(SRID_NORTH_LAMBERT);
	}

	/* Are these data antarctic? Lambert Azimuthal Equal Area South. */
	if ( center.y < -70.0 && ywidth < 45.0 )
	{
		PG_RETURN_INT32(SRID_SOUTH_LAMBERT);
	}

	/*
	** Can we fit these data into one UTM zone?
	** We will assume we can push things as
	** far as a half zone past a zone boundary.
	** Note we have no handling for the date line in here.
	*/
	if ( xwidth < 6.0 )
	{
		int zone = floor((center.x + 180.0) / 6.0);

		if ( zone > 59 ) zone = 59;

		/* Are these data below the equator? UTM South. */
		if ( center.y < 0.0 )
		{
			PG_RETURN_INT32( SRID_SOUTH_UTM_START + zone );
		}
		/* Are these data above the equator? UTM North. */
		else
		{
			PG_RETURN_INT32( SRID_NORTH_UTM_START + zone );
		}
	}

	/*
	** Can we fit into a custom LAEA area? (30 degrees high, variable width)
	** We will allow overlap into adjoining areas, but use a slightly narrower test (25) to try
	** and minimize the worst case.
	** Again, we are hoping the dateline doesn't trip us up much
	*/
	if ( ywidth < 25.0 )
	{
		int xzone = -1;
		int yzone = 3 + floor(center.y / 30.0); /* (range of 0-5) */

		/* Equatorial band, 12 zones, 30 degrees wide */
		if ( (yzone == 2 || yzone == 3) && xwidth < 30.0 )
		{
			xzone = 6 + floor(center.x / 30.0);
		}
		/* Temperate band, 8 zones, 45 degrees wide */
		else if ( (yzone == 1 || yzone == 4) && xwidth < 45.0 )
		{
			xzone = 4 + floor(center.x / 45.0);
		}
		/* Arctic band, 4 zones, 90 degrees wide */
		else if ( (yzone == 0 || yzone == 5) && xwidth < 90.0 )
		{
			xzone = 2 + floor(center.x / 90.0);
		}

		/* Did we fit into an appropriate xzone? */
		if ( xzone != -1 )
		{
			PG_RETURN_INT32(SRID_LAEA_START + 20 * yzone + xzone);
		}
	}

	/*
	** Running out of options... fall-back to Mercator
	** and hope for the best.
	*/
	PG_RETURN_INT32(SRID_WORLD_MERCATOR);

}

/*
** geography_project(GSERIALIZED *g, distance, azimuth)
** returns point of projection given start point,
** azimuth in radians (bearing) and distance in meters
*/
PG_FUNCTION_INFO_V1(geography_project);
Datum geography_project(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom = NULL;
	LWPOINT *lwp_projected;
	GSERIALIZED *g = NULL;
	GSERIALIZED *g_out = NULL;
	double azimuth = 0.0;
	double distance;
	SPHEROID s;
	uint32_t type;

	/* Return NULL on NULL distance or geography */
	if ( PG_NARGS() < 2 || PG_ARGISNULL(0) || PG_ARGISNULL(1) )
		PG_RETURN_NULL();

	/* Get our geometry object loaded into memory. */
	g = PG_GETARG_GSERIALIZED_P(0);

	/* Only return for points. */
	type = gserialized_get_type(g);
	if ( type != POINTTYPE )
	{
		elog(ERROR, "ST_Project(geography) is only valid for point inputs");
		PG_RETURN_NULL();
	}

	distance = PG_GETARG_FLOAT8(1); /* Distance in Meters */
	lwgeom = lwgeom_from_gserialized(g);

	/* EMPTY things cannot be projected from */
	if ( lwgeom_is_empty(lwgeom) )
	{
		lwgeom_free(lwgeom);
		elog(ERROR, "ST_Project(geography) cannot project from an empty start point");
		PG_RETURN_NULL();
	}

	if ( PG_NARGS() > 2 && ! PG_ARGISNULL(2) )
		azimuth = PG_GETARG_FLOAT8(2); /* Azimuth in Radians */

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g), &s);

	/* Handle the zero distance case */
	if( FP_EQUALS(distance, 0.0) )
	{
		PG_RETURN_POINTER(g);
	}

	/* Calculate the length */
	lwp_projected = lwgeom_project_spheroid(lwgeom_as_lwpoint(lwgeom), &s, distance, azimuth);

	/* Something went wrong... */
	if ( lwp_projected == NULL )
	{
		elog(ERROR, "lwgeom_project_spheroid returned null");
		PG_RETURN_NULL();
	}

	/* Clean up, but not all the way to the point arrays */
	lwgeom_free(lwgeom);
	g_out = geography_serialize(lwpoint_as_lwgeom(lwp_projected));
	lwpoint_free(lwp_projected);

	PG_FREE_IF_COPY(g, 0);
	PG_RETURN_POINTER(g_out);
}


/*
** geography_azimuth(GSERIALIZED *g1, GSERIALIZED *g2)
** returns direction between points (north = 0)
** azimuth (bearing) and distance
*/
PG_FUNCTION_INFO_V1(geography_azimuth);
Datum geography_azimuth(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom1 = NULL;
	LWGEOM *lwgeom2 = NULL;
	GSERIALIZED *g1 = NULL;
	GSERIALIZED *g2 = NULL;
	double azimuth;
	SPHEROID s;
	uint32_t type1, type2;

	/* Get our geometry object loaded into memory. */
	g1 = PG_GETARG_GSERIALIZED_P(0);
	g2 = PG_GETARG_GSERIALIZED_P(1);

	/* Only return for points. */
	type1 = gserialized_get_type(g1);
	type2 = gserialized_get_type(g2);
	if ( type1 != POINTTYPE || type2 != POINTTYPE )
	{
		elog(ERROR, "ST_Azimuth(geography, geography) is only valid for point inputs");
		PG_RETURN_NULL();
	}

	lwgeom1 = lwgeom_from_gserialized(g1);
	lwgeom2 = lwgeom_from_gserialized(g2);

	/* EMPTY things cannot be used */
	if ( lwgeom_is_empty(lwgeom1) || lwgeom_is_empty(lwgeom2) )
	{
		lwgeom_free(lwgeom1);
		lwgeom_free(lwgeom2);
		elog(ERROR, "ST_Azimuth(geography, geography) cannot work with empty points");
		PG_RETURN_NULL();
	}

	/* Initialize spheroid */
	spheroid_init_from_srid(fcinfo, gserialized_get_srid(g1), &s);

	/* Calculate the direction */
	azimuth = lwgeom_azumith_spheroid(lwgeom_as_lwpoint(lwgeom1), lwgeom_as_lwpoint(lwgeom2), &s);

	/* Clean up */
	lwgeom_free(lwgeom1);
	lwgeom_free(lwgeom2);

	PG_FREE_IF_COPY(g1, 0);
	PG_FREE_IF_COPY(g2, 1);

	/* Return NULL for unknown (same point) azimuth */
	if( isnan(azimuth) )
	{
		PG_RETURN_NULL();
	}

	PG_RETURN_FLOAT8(azimuth);
}



/*
** geography_segmentize(GSERIALIZED *g1, double max_seg_length)
** returns densified geometry with no segment longer than max
*/
PG_FUNCTION_INFO_V1(geography_segmentize);
Datum geography_segmentize(PG_FUNCTION_ARGS)
{
	LWGEOM *lwgeom1 = NULL;
	LWGEOM *lwgeom2 = NULL;
	GSERIALIZED *g1 = NULL;
	GSERIALIZED *g2 = NULL;
	double max_seg_length;
	uint32_t type1;

	/* Get our geometry object loaded into memory. */
	g1 = PG_GETARG_GSERIALIZED_P(0);
	type1 = gserialized_get_type(g1);

	/* Convert max_seg_length from metric units to radians */
	max_seg_length = PG_GETARG_FLOAT8(1) / WGS84_RADIUS;

	/* We can't densify points or points, reflect them back */
	if ( type1 == POINTTYPE || type1 == MULTIPOINTTYPE || gserialized_is_empty(g1) )
		PG_RETURN_POINTER(g1);

	/* Deserialize */
	lwgeom1 = lwgeom_from_gserialized(g1);

	/* Calculate the densified geometry */
	lwgeom2 = lwgeom_segmentize_sphere(lwgeom1, max_seg_length);

	/*
	** Set the geodetic flag so subsequent
	** functions do the right thing.
	*/
	lwgeom_set_geodetic(lwgeom2, true);

	/* Recalculate the boxes after re-setting the geodetic bit */
	lwgeom_drop_bbox(lwgeom2);

	/* We are trusting geography_serialize will add a box if needed */
	g2 = geography_serialize(lwgeom2);

	/* Clean up */
	lwgeom_free(lwgeom1);
	lwgeom_free(lwgeom2);
	PG_FREE_IF_COPY(g1, 0);

	PG_RETURN_POINTER(g2);
}


/* Interpolate a point along a geographic line. */
void
interpolate_point4d_sphere(
        const POINT3D *p1, const POINT3D *p2,
        const POINT4D *v1, const POINT4D *v2,
        double f,
        POINT4D *p)
{
    /* Calculate interpolated point */
    POINT3D mid;
    mid.x = p1->x + ((p2->x - p1->x) * f);
    mid.y = p1->y + ((p2->y - p1->y) * f);
    mid.z = p1->z + ((p2->z - p1->z) * f);
    normalize(&mid);

    /* Calculate z/m values */
    GEOGRAPHIC_POINT g;
    cart2geog(&mid, &g);
    p->x = rad2deg(g.lon);
    p->y = rad2deg(g.lat);
    p->z = v1->z + ((v2->z - v1->z) * f);
    p->m = v1->m + ((v2->m - v1->m) * f);
}

double ptarray_length_sphere(const POINTARRAY *pa)
{
    GEOGRAPHIC_POINT a, b;
    POINT4D p;
    uint32_t i;
    double length = 0.0;

    /* Return zero on non-sensical inputs */
    if ( ! pa || pa->npoints < 2 )
        return 0.0;

    /* Initialize first point */
    getPoint4d_p(pa, 0, &p);
    geographic_point_init(p.x, p.y, &a);

    /* Loop and sum the length for each segment */
    for ( i = 1; i < pa->npoints; i++ )
    {
        getPoint4d_p(pa, i, &p);
        geographic_point_init(p.x, p.y, &b);
        /* Add this segment length to the total */
        length +=  sphere_distance(&a, &b);
    }
    return length;
}

POINTARRAY* lwline_interpolate_points_sphere(const LWLINE *line, double length_fraction,
                                             char repeat)
{
    POINT4D pt;
    uint32_t i;
    uint32_t points_to_interpolate;
    uint32_t points_found = 0;
    double length;
    double length_fraction_increment = length_fraction;
    double length_fraction_consumed = 0;
    char has_z = (char) lwgeom_has_z(lwline_as_lwgeom(line));
    char has_m = (char) lwgeom_has_m(lwline_as_lwgeom(line));
    const POINTARRAY* ipa = line->points;
    POINTARRAY* opa;
    POINT4D p1, p2;
    POINT3D q1, q2;
    GEOGRAPHIC_POINT g1, g2;

    /* Empty.InterpolatePoint == Point Empty */
    if ( lwline_is_empty(line) )
    {
        return ptarray_construct_empty(has_z, has_m, 0);
    }

    /* If distance is one of the two extremes, return the point on that
     * end rather than doing any computations
     */
    if ( length_fraction == 0.0 || length_fraction == 1.0 )
    {
        if ( length_fraction == 0.0 )
            getPoint4d_p(ipa, 0, &pt);
        else
            getPoint4d_p(ipa, ipa->npoints-1, &pt);

        opa = ptarray_construct(has_z, has_m, 1);
        ptarray_set_point4d(opa, 0, &pt);

        return opa;
    }

    /* Interpolate points along the line */
    length = ptarray_length_sphere(ipa);
    points_to_interpolate = repeat ? (uint32_t) floor(1 / length_fraction) : 1;
    opa = ptarray_construct(has_z, has_m, points_to_interpolate);

    getPoint4d_p(ipa, 0, &p1);
    geographic_point_init(p1.x, p1.y, &g1);
    for ( i = 0; i < ipa->npoints - 1 && points_found < points_to_interpolate; i++ )
    {
        getPoint4d_p(ipa, i+1, &p2);
        geographic_point_init(p2.x, p2.y, &g2);
        double segment_length_frac = sphere_distance(&g1, &g2) / length;

        /* If our target distance is before the total length we've seen
         * so far. create a new point some distance down the current
         * segment.
         */
        while ( length_fraction < length_fraction_consumed + segment_length_frac && points_found < points_to_interpolate )
        {
            double segment_fraction = (length_fraction - length_fraction_consumed) / segment_length_frac;
            geog2cart(&g1, &q1);
            geog2cart(&g2, &q2);
            interpolate_point4d_sphere(&q1, &q2, &p1, &p2, segment_fraction, &pt);
            ptarray_set_point4d(opa, points_found++, &pt);
            length_fraction += length_fraction_increment;
        }

        length_fraction_consumed += segment_length_frac;

        p1 = p2;
        g1 = g2;
    }

    /* Return the last point on the line. This shouldn't happen, but
     * could if there's some floating point rounding errors. */
    if (points_found < points_to_interpolate) {
        getPoint4d_p(ipa, ipa->npoints - 1, &pt);
        ptarray_set_point4d(opa, points_found, &pt);
    }

    return opa;
}

/*
** geography_line_interpolate_point(GSERIALIZED *g1, double distance_fraction, bool repeat)
** returns a point interpolated along a geographic line
*/
PG_FUNCTION_INFO_V1(geography_line_interpolate_point);
Datum geography_line_interpolate_point(PG_FUNCTION_ARGS)
{
    GSERIALIZED *g1 = PG_GETARG_GSERIALIZED_P(0);
    double distance_fraction = PG_GETARG_FLOAT8(1);
    bool repeat = PG_NARGS() > 2 && PG_GETARG_BOOL(2);
    int srid = gserialized_get_srid(g1);
    LWLINE* lwline;
    LWGEOM* lwresult;
    POINTARRAY* opa;
    GSERIALIZED *result;

    if ( distance_fraction < 0 || distance_fraction > 1 )
    {
        elog(ERROR,"line_interpolate_point: 2nd arg isn't within [0,1]");
        PG_FREE_IF_COPY(g1, 0);
        PG_RETURN_NULL();
    }

    if ( gserialized_get_type(g1) != LINETYPE )
    {
        elog(ERROR,"line_interpolate_point: 1st arg isn't a line");
        PG_FREE_IF_COPY(g1, 0);
        PG_RETURN_NULL();
    }

    lwline = lwgeom_as_lwline(lwgeom_from_gserialized(g1));
    opa = lwline_interpolate_points_sphere(lwline, distance_fraction, repeat);

    lwgeom_free(lwline_as_lwgeom(lwline));
    PG_FREE_IF_COPY(g1, 0);

    if (opa->npoints <= 1)
    {
        lwresult = lwpoint_as_lwgeom(lwpoint_construct(srid, NULL, opa));
    } else {
        lwresult = lwmpoint_as_lwgeom(lwmpoint_construct(srid, opa));
    }

    lwgeom_set_geodetic(lwresult, true);
    result = geography_serialize(lwresult);
    lwgeom_free(lwresult);

    PG_RETURN_POINTER(result);
}


/* Locate a point along a geographic line. */
double
ptarray_locate_point_sphere(const POINTARRAY *pa, const POINT4D *p,
                            double tolerance, double *mindistout, POINT4D *closest)
{
    GEOGRAPHIC_EDGE e;
    GEOGRAPHIC_POINT a, b, nearest;
    POINT4D p1, p2;
    POINT2D p2d;
    uint32_t i, seg = 0;
    double distance, result;
    long double fraction, /* Used for computing Z and M values of the closest point */
            length, /* Length of the current segment */
            seglength, /* length of the segment where the closest point is located */
            partlength = 0.0, /* length from the beginning of the point array to the closest point */
            totlength = 0.0;  /* length of the point array */

    /* Initialize target point */
    geographic_point_init(p->x, p->y, &a);

    /* Handle point/point case here */
    if ( pa->npoints <= 1)
    {
        if ( pa->npoints == 1 && mindistout )
        {
            getPoint4d_p(pa, 0, &p1);
            geographic_point_init(p1.x, p1.y, &b);
            *mindistout = sphere_distance(&a, &b);
        }
        return 0.0;
    }

    /* Make distance really big, so that everything will be smaller than it */
    distance = FLT_MAX;

    /* Initialize first point of array */
    getPoint4d_p(pa, 0, &p1);
    geographic_point_init(p1.x, p1.y, &(e.start));

    /* Iterate through the edges in the line */
    for ( i = 1; i < pa->npoints; i++ )
    {
        double d;
        getPoint4d_p(pa, i, &p2);
        geographic_point_init(p2.x, p2.y, &(e.end));
        /* Get the spherical distance between point and edge */
        d = edge_distance_to_point(&e, &a, &b);
        /* New shortest distance! Record this distance/location/segment */
        if ( d < distance )
        {
            distance = d;
            nearest = b;
            seg = i - 1;
        }
        /* We've gotten closer than the tolerance... */
        if ( d < tolerance )
            break;

        e.start = e.end;
    }

    /* Initialize first point of array */
    getPoint4d_p(pa, 0, &p1);
    geographic_point_init(p1.x, p1.y, &a);

    /* Loop and sum the length for each segment */
    for ( i = 1; i < pa->npoints; i++ )
    {
        getPoint4d_p(pa, i, &p2);
        geographic_point_init(p2.x, p2.y, &b);

        /* Compute length of current segment */
        length = sphere_distance(&a, &b);

        /* Add segment length to the partial and total length */
        if (i < seg)
            partlength += length;
        else if (i == seg)
            seglength = length;
        totlength += length;

        /* B gets incremented in the next loop, so we save the value here */
        a = b;
    }

    /* Get the points defining the segment of the closest point */
    getPoint4d_p(pa, seg, &p1);
    getPoint4d_p(pa, seg + 1, &p2);

    /* Compute distance from beginning of the segment to closest point */
    geographic_point_init(p1.x, p1.y, &a);
    length = sphere_distance(&a, &nearest);

    /* Add this length to the partial length */
    partlength += length;

    /* Set output parameters */
    if ( mindistout )
        *mindistout = distance;
    if ( closest )
    {
        /* Set lon and lat for output parameter */
        closest->x = p2d.x = rad2deg(nearest.lon);
        closest->y = p2d.y = rad2deg(nearest.lat);

        /* For robustness, return the original point in line when
         * closest point ~= one of the points in line */
        if (p2d_same(&p2d, getPoint2d_cp(pa, seg)))
            getPoint4d_p(pa, seg, closest);
        else if (p2d_same(&p2d, getPoint2d_cp(pa, seg + 1)))
            getPoint4d_p(pa, seg + 1, closest);
        else
        {
            if (ptarray_has_z(pa) || ptarray_has_m(pa))
            {
                fraction = length / seglength;
                closest->z = p1.z + (double) ((long double) (p2.z - p1.z) * fraction);
                closest->m = p1.m + (double) ((long double) (p2.m - p1.m) * fraction);
            }
            else
            {
                closest->z = NO_Z_VALUE;
                closest->m = NO_M_VALUE;
            }
        }
    }

    /* Location of any point on a zero-length line is 0 */
    /* See http://trac.osgeo.org/postgis/ticket/1772#comment:2 */
    if ( totlength == 0 )
        return 0.0;

    /* For robustness, force 0/1 when closest point == start/endpoint */
    getPoint4d_p(pa, 0, &p1);
    getPoint4d_p(pa, pa->npoints - 1, &p2);
    if ( seg == 0 && p4d_same(closest, &p1) )
        return 0.0;
    if ( seg >= (pa->npoints - 2) && p4d_same(closest, &p2) )
        return 1.0;

    result = (double) (partlength / totlength);
    return result;
}

/*
** geography_line_locate_point(GSERIALIZED *g1, GSERIALIZED *g2)
** returns a float between 0 and 1 representing the location of the closest point on LineString to the given Point
*/
PG_FUNCTION_INFO_V1(geography_line_locate_point);
Datum geography_line_locate_point(PG_FUNCTION_ARGS)
{
    GSERIALIZED *g1 = PG_GETARG_GSERIALIZED_P(0);
    GSERIALIZED *g2 = PG_GETARG_GSERIALIZED_P(1);
    LWLINE *lwline;
    LWPOINT *lwpoint;
    POINT4D p, proj;
    double ret;

    if ( gserialized_get_type(g1) != LINETYPE )
    {
        elog(ERROR,"line_locate_point: 1st arg isn't a line");
        PG_RETURN_NULL();
    }
    if ( gserialized_get_type(g2) != POINTTYPE )
    {
        elog(ERROR,"line_locate_point: 2st arg isn't a point");
        PG_RETURN_NULL();
    }

    gserialized_error_if_srid_mismatch(g1, g2, __func__);

    lwline = lwgeom_as_lwline(lwgeom_from_gserialized(g1));
    lwpoint = lwgeom_as_lwpoint(lwgeom_from_gserialized(g2));

    lwpoint_getPoint4d_p(lwpoint, &p);

    ret = ptarray_locate_point_sphere(lwline->points, &p, FP_TOLERANCE, NULL, &proj);

    PG_RETURN_FLOAT8(ret);
}



/*
** geography_closestpoint(GSERIALIZED *g1, GSERIALIZED *g2)
** returns the point in first input geography that is closest to the second input geography in 2d
*/
PG_FUNCTION_INFO_V1(geography_closestpoint);
Datum geography_closestpoint(PG_FUNCTION_ARGS)
{
    GSERIALIZED* g1 = NULL;
    GSERIALIZED* g2 = NULL;
    LWGEOM *point;
    GSERIALIZED* result;

    /* Get our geography objects loaded into memory. */
    g1 = PG_GETARG_GSERIALIZED_P(0);
    g2 = PG_GETARG_GSERIALIZED_P(1);

    gserialized_error_if_srid_mismatch(g1, g2, __func__);

    /* Return NULL on empty arguments. */
    if ( gserialized_is_empty(g1) || gserialized_is_empty(g2) )
    {
        PG_FREE_IF_COPY(g1, 0);
        PG_FREE_IF_COPY(g2, 1);
        PG_RETURN_NULL();
    }

    point = geography_tree_closestpoint(g1, g2, FP_TOLERANCE);

    if (lwgeom_is_empty(point))
        PG_RETURN_NULL();

    result = geography_serialize(point);
    lwgeom_free(point);

    PG_FREE_IF_COPY(g1, 0);
    PG_FREE_IF_COPY(g2, 1);
    PG_RETURN_POINTER(result);
}



/*
** geography_shortestline(GSERIALIZED *g1, GSERIALIZED *g2)
** returns the 2-dimensional shortest line between two geographies
*/
PG_FUNCTION_INFO_V1(geography_shortestline);
Datum geography_shortestline(PG_FUNCTION_ARGS)
{
    GSERIALIZED* g1 = NULL;
    GSERIALIZED* g2 = NULL;
    LWGEOM *line;
    GSERIALIZED* result;

    /* Get our geography objects loaded into memory. */
    g1 = PG_GETARG_GSERIALIZED_P(0);
    g2 = PG_GETARG_GSERIALIZED_P(1);

    gserialized_error_if_srid_mismatch(g1, g2, __func__);

    /* Return NULL on empty arguments. */
    if ( gserialized_is_empty(g1) || gserialized_is_empty(g2) )
    {
        PG_FREE_IF_COPY(g1, 0);
        PG_FREE_IF_COPY(g2, 1);
        PG_RETURN_NULL();
    }

    line = geography_tree_shortestline(g1, g2, FP_TOLERANCE);

    if (lwgeom_is_empty(line))
        PG_RETURN_NULL();

    result = geography_serialize(line);
    lwgeom_free(line);

    PG_FREE_IF_COPY(g1, 0);
    PG_FREE_IF_COPY(g2, 1);
    PG_RETURN_POINTER(result);
}


/*
** New version of ST_Segmentize that allows to perform the calculation on a sphere or spheroid
** @param p1, p2 - 3-space points we are interpolating between
** @param v1, v2 - real values and z/m values
** @param d, max_seg_length - current segment length and segment limit
*/
static int ptarray_segmentize_spheroid_edge_recursive (
        const POINT3D *p1, const POINT3D *p2,
        const POINT4D *v1, const POINT4D *v2,
        double d, double max_seg_length,
        POINTARRAY *pa)
{
    GEOGRAPHIC_POINT g;
    /* Reached the terminal leaf in recursion. Add */
    /* the left-most point to the pointarray here */
    /* We recurse down the left side first, so outputs should */
    /* end up added to the array in order this way */
    if (d <= max_seg_length)
    {
        POINT4D p;
        cart2geog(p1, &g);
        p.x = v1->x;
        p.y = v1->y;
        p.z = v1->z;
        p.m = v1->m;
        return ptarray_append_point(pa, &p, LW_FALSE);
    }
        /* Find the mid-point and recurse on the left and then the right */
    else
    {
        /* Calculate mid-point */
        POINT3D mid;
        mid.x = (p1->x + p2->x) / 2.0;
        mid.y = (p1->y + p2->y) / 2.0;
        mid.z = (p1->z + p2->z) / 2.0;
        normalize(&mid);

        /* Calculate z/m mid-values */
        POINT4D midv;
        cart2geog(&mid, &g);
        midv.x = rad2deg(g.lon);
        midv.y = rad2deg(g.lat);
        midv.z = (v1->z + v2->z) / 2.0;
        midv.m = (v1->m + v2->m) / 2.0;
        /* Recurse on the left first */
        ptarray_segmentize_spheroid_edge_recursive(p1, &mid, v1, &midv, d/2.0, max_seg_length, pa);
        ptarray_segmentize_spheroid_edge_recursive(&mid, p2, &midv, v2, d/2.0, max_seg_length, pa);
        return LW_SUCCESS;
    }
}

/*
** Create a new point array with no segment longer than the input segment length (expressed in radians!)
** @param pa_in - input point array pointer
** @param max_seg_length - maximum output segment length in radians
*/
static POINTARRAY*
ptarray_segmentize_spheroid(const POINTARRAY *pa_in, double max_seg_length, const SPHEROID *s)
{
    POINTARRAY *pa_out;
    int hasz = ptarray_has_z(pa_in);
    int hasm = ptarray_has_m(pa_in);
    POINT4D p1, p2;
    POINT3D q1, q2;
    GEOGRAPHIC_POINT g1, g2;
    uint32_t i;

    /* Just crap out on crazy input */
    if ( ! pa_in )
        lwerror("%s: null input pointarray", __func__);
    if ( max_seg_length <= 0.0 )
        lwerror("%s: maximum segment length must be positive", __func__);

    /* Empty starting array */
    pa_out = ptarray_construct_empty(hasz, hasm, pa_in->npoints);

    /* Simple loop per edge */
    for (i = 1; i < pa_in->npoints; i++)
    {
        getPoint4d_p(pa_in, i-1, &p1);
        getPoint4d_p(pa_in, i, &p2);
        geographic_point_init(p1.x, p1.y, &g1);
        geographic_point_init(p2.x, p2.y, &g2);

        /* Skip duplicate points (except in case of 2-point lines!) */
        if ((pa_in->npoints > 2) && p4d_same(&p1, &p2))
            continue;

        /* How long is this edge? */
        double d;
        /* Special sphere case */
        if ( s->a == s->b )
            d = s->radius * sphere_distance(&g1, &g2);
            /* Spheroid case */
        else
            d = spheroid_distance(&g1, &g2, s);

        if (d > max_seg_length)
        {
            geog2cart(&g1, &q1);
            geog2cart(&g2, &q2);
            /* 3-d end points, XYZM end point, current edge size, min edge size */
            ptarray_segmentize_spheroid_edge_recursive(&q1, &q2, &p1, &p2, d, max_seg_length, pa_out);
        }
            /* If we don't segmentize, we need to add first point manually */
        else
        {
            ptarray_append_point(pa_out, &p1, LW_TRUE);
        }
    }
    /* Always add the last point */
    ptarray_append_point(pa_out, &p2, LW_TRUE);
    return pa_out;
}

/*
** Create a new, densified geometry where no segment is longer than max_seg_length.
** Input geometry is not altered, output geometry must be freed by caller.
** @param lwg_in = input geometry
** @param max_seg_length = maximum segment length in radians
*/
LWGEOM*
lwgeom_segmentize_spheroid(const LWGEOM *lwg_in, double max_seg_length, const SPHEROID *s)
{
    POINTARRAY *pa_out;
    LWLINE *lwline;
    LWPOLY *lwpoly_in, *lwpoly_out;
    LWCOLLECTION *lwcol_in, *lwcol_out;
    uint32_t i;

    /* Reflect NULL */
    if ( ! lwg_in )
        return NULL;

    /* Clone empty */
    if ( lwgeom_is_empty(lwg_in) )
        return lwgeom_clone(lwg_in);

    switch (lwg_in->type)
    {
        case MULTIPOINTTYPE:
        case POINTTYPE:
            return lwgeom_clone_deep(lwg_in);
            break;
        case LINETYPE:
            lwline = lwgeom_as_lwline(lwg_in);
            pa_out = ptarray_segmentize_spheroid(lwline->points, max_seg_length, s);
            return lwline_as_lwgeom(lwline_construct(lwg_in->srid, NULL, pa_out));
            break;
        case POLYGONTYPE:
            lwpoly_in = lwgeom_as_lwpoly(lwg_in);
            lwpoly_out = lwpoly_construct_empty(lwg_in->srid, lwgeom_has_z(lwg_in), lwgeom_has_m(lwg_in));
            for ( i = 0; i < lwpoly_in->nrings; i++ )
            {
                pa_out = ptarray_segmentize_spheroid(lwpoly_in->rings[i], max_seg_length, s);
                lwpoly_add_ring(lwpoly_out, pa_out);
            }
            return lwpoly_as_lwgeom(lwpoly_out);
            break;
        case MULTILINETYPE:
        case MULTIPOLYGONTYPE:
        case COLLECTIONTYPE:
            lwcol_in = lwgeom_as_lwcollection(lwg_in);
            lwcol_out = lwcollection_construct_empty(lwg_in->type, lwg_in->srid, lwgeom_has_z(lwg_in), lwgeom_has_m(lwg_in));
            for ( i = 0; i < lwcol_in->ngeoms; i++ )
            {
                lwcollection_add_lwgeom(lwcol_out, lwgeom_segmentize_spheroid(lwcol_in->geoms[i], max_seg_length, s));
            }
            return lwcollection_as_lwgeom(lwcol_out);
            break;
        default:
            lwerror("lwgeom_segmentize_spheroid: unsupported input geometry type: %d - %s",
                    lwg_in->type, lwtype_name(lwg_in->type));
            break;
    }

    lwerror("lwgeom_segmentize_spheroid got to the end of the function, should not happen");
    return NULL;
}

/*
** geography_segmentize1(GSERIALIZED *g1, double max_seg_length, boolean use_spheroid)
** returns densified geometry with no segment longer than max
*/
PG_FUNCTION_INFO_V1(geography_segmentize1);
Datum geography_segmentize1(PG_FUNCTION_ARGS)
{
    LWGEOM *lwgeom1 = NULL;
    LWGEOM *lwgeom2 = NULL;
    GSERIALIZED *g1 = NULL;
    GSERIALIZED *g2 = NULL;
    double max_seg_length;
    bool use_spheroid = true;
    uint32_t type1;
    SPHEROID s;

    /* Get our geometry object loaded into memory. */
    g1 = PG_GETARG_GSERIALIZED_P(0);
    type1 = gserialized_get_type(g1);

    /* Get max_seg_length in metric units */
    max_seg_length = PG_GETARG_FLOAT8(1);

    /* Read calculation type */
    if ( PG_NARGS() > 2 && ! PG_ARGISNULL(2) )
        use_spheroid = PG_GETARG_BOOL(2);

    /* We can't densify points or points, reflect them back */
    if ( type1 == POINTTYPE || type1 == MULTIPOINTTYPE || gserialized_is_empty(g1) )
        PG_RETURN_POINTER(g1);

    /* Initialize spheroid */
    /* We currently cannot use the following statement since PROJ4 API is not
     * available directly to MobilityDB. */
    // spheroid_init_from_srid(fcinfo, srid, &s);
    spheroid_init(&s, WGS84_MAJOR_AXIS, WGS84_MINOR_AXIS);

    /* Set to sphere if requested */
    if ( ! use_spheroid )
        s.a = s.b = s.radius;

    /* Deserialize */
    lwgeom1 = lwgeom_from_gserialized(g1);

    /* Calculate the densified geometry */
    lwgeom2 = lwgeom_segmentize_spheroid(lwgeom1, max_seg_length, &s);

    /*
    ** Set the geodetic flag so subsequent
    ** functions do the right thing.
    */
    lwgeom_set_geodetic(lwgeom2, true);

    /* Recalculate the boxes after re-setting the geodetic bit */
    lwgeom_drop_bbox(lwgeom2);

    /* We are trusting geography_serialize will add a box if needed */
    g2 = geography_serialize(lwgeom2);

    /* Clean up */
    lwgeom_free(lwgeom1);
    lwgeom_free(lwgeom2);
    PG_FREE_IF_COPY(g1, 0);

    PG_RETURN_POINTER(g2);
}
