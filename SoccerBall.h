/* Soccerball, Copyright (c) 2004-2020 Yuriy Yakimenko
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation.  No representations are made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or 
 * implied warranty.
 */

#pragma once

#include "math.h"
#include "assert.h"

#include <vector>
#include <X11/Xlib.h>

#include "GenGraphics.h"

#define BACKGD  0xA0A0C0  // background
#define DKGRAY  0x606090  // inner part of pentagons 

#define BLACK 	0x000000
#define WHITE 	0xFFFFFF

#define BALL_RADIUS 190
#define SEAM_THICKNESS 0.01

class DPoint
{
	public:
	double x,y,z;

	public:

	DPoint () : x(0), y(0), z(0) { }

	DPoint (double _x, double _y, double _z)
	{
		x = _x; 
		y = _y;
		z = _z;
	}
	
	void Copy (const DPoint & point)
	{
		x = point.x;
		y = point.y;
		z = point.z;
	}  

	void Assign (double _x, double _y, double _z)
	{
		x = _x; 
		y = _y;
		z = _z;
	}

	void Normalize ()
	{
		double k = 1.0 / sqrt(x*x + y*y + z*z);
		x = x * k;
		y = y * k;
		z = z * k;
	}

	double SquaredDistance (const DPoint & other) const 
	{
		return	(x - other.x) * (x - other.x) +
				(y - other.y) * (y - other.y) +
				(z - other.z) * (z - other.z) ;
	}
};


class GeneralPath
{
	private:
	std::vector<std::pair<double,double> > path;
	public:
	GeneralPath() {}
	~GeneralPath()
	{
		path.clear();
	}

	void lineTo (double x, double y)
	{
		path.push_back(std::make_pair(x,y));
	}

	void moveTo (double x, double y)
	{
		path.clear();
		path.push_back(std::make_pair(x,y));
	}

	void closePath() {}

	void fill (GC & gc, GeneralGraphics * gg) const 
	{
		gg->fillPolygon (gc, path);
	}

};

class SoccerBall 
{
	private:

	GeneralGraphics * gg;

	float radius;
	int centerX, centerY;
	int iWidth, iHeight;
	DPoint * stored;
	DPoint * vertex;
	int seams_parts;
	double delta, tilt, delta1, tilt1;
	bool captured;

	public :

	bool mouseDown (int x, int y)
	{
		double x1 = x - centerX;
		double y1 = y - centerY;

		if (sqrt(x1*x1 + y1*y1) < radius)
		{
			delta1 = asin (x1/sqrt(radius*radius-y1*y1));
			tilt1 = asin (y1/sqrt(radius*radius-x1*x1));

			captured = true;
		}

		return true;
	}

	bool mouseUp (int x, int y)
	{
		if (captured)
		{
			captured = false;
			updateStorage (delta, tilt);
			delta = 0;
			tilt = 0;
		}

		return true;
	}

	bool mouseDrag (int x, int y)
	{
		double x1 = x - centerX;
		double y1 = y - centerY;

		if (captured && sqrt(x1*x1 + y1*y1) < radius)
		{
			double delta2 = asin (x1/sqrt (radius*radius-y1*y1));
			double  tilt2 = asin (y1/sqrt(radius*radius-x1*x1));

			delta = (delta2 - delta1); 
			tilt = (tilt2 - tilt1); 

			updateCurrent (delta, tilt);
		}

		return true;
	}

	private:

	void TransformDelta (DPoint & point, double delta)
	{
		double x = point.x;
		double z = (point.z) ;

		point.x = x * cos(delta) + z * sin(delta) ;
		point.z = z * cos(delta) - x * sin(delta) ;
	}

	void TransformTilt (DPoint & point, double tilt)
	{
		double y = point.y ;
		double z = point.z ;

		point.y = y * cos (tilt) + z * sin(tilt) ;
		point.z = z * cos (tilt) - y * sin(tilt) ;
	}

	double GetZCoord (const DPoint & pt)
	{
		return pt.z;
	}

	double SquaredDistance (const DPoint & pt1, const DPoint & pt2)
	{
		return pt1.SquaredDistance (pt2);
	}

	double AngleDPoints (const DPoint & point1, const DPoint & point2)
	{
		double d = (point1.x * point2.x)  +
				 	(point1.y * point2.y)  +
					(point1.z * point2.z) ;

		return asin(d);
	}

	void updateStorage (double delta, double tilt)
	{
		for (int i = 0; i < 12; i++)
		{
			TransformDelta (stored[i], delta);
			TransformTilt (stored[i], tilt);
		}
	}

	void updateCurrent (double delta, double tilt)
	{
		for (int i = 0; i < 12; i++)
		{
			vertex[i].Copy(stored[i]);
			TransformDelta (vertex[i], delta);
			TransformTilt (vertex[i], tilt);
		}
	}

	void initVertices ()
	{
		double alpha;
		double beta = atan (2.0);

		int i;

		vertex [0].x = 0;
		vertex [0].y = 1;
		vertex [0].z = 0;

		for (i = 1; i <= 5; i++)
		{
			alpha = 0.4 * M_PI * (i-1); 
			vertex[i].x = sin(beta) * cos(alpha);
			vertex[i].y = cos(beta);
			vertex[i].z = sin(beta) * sin(alpha);
		}

		for (i = 6; i < 12; i++)
		{
			vertex[i].x = -vertex[i-6].x;
			vertex[i].y = -vertex[i-6].y;
			vertex[i].z = -vertex[i-6].z;
		}	

		for (i = 0; i < 12; i++)
		{
			stored[i].Copy(vertex[i]);
		}
	}

	void FindIntermediatePoint (const DPoint & pt1, const DPoint & pt2, DPoint & out, double fraction)
	{
		out.x = pt1.x + (pt2.x - pt1.x) * fraction;
		out.y = pt1.y + (pt2.y - pt1.y) * fraction;
		out.z = pt1.z + (pt2.z - pt1.z) * fraction;

		out.Normalize();
	}

	int findClosest (const DPoint & from, int start, const DPoint * points, int cnt, double min)
	{
		double dist = 0;

		for (int i = start; i < cnt; i++)
			{
			dist = from.SquaredDistance (points[i]);

			if (dist < min)	return i;
			}
		return -1;
	}

	void drawPaintedPentagons (const DPoint * pentagon_points, GC & gc)
	{
		DPoint * pent_set1 = new DPoint [5];
		int j, cnt, index, closest, total_parts = 5 * seams_parts;
		DPoint temp;
		DPoint * polyline = new DPoint[total_parts];
		DPoint * edge = new DPoint[2];

		for (int i = 0; i < 12; i++)
		{
			closest = 0;
			cnt=0;
			while (0 <= (closest = findClosest (vertex[i], closest, pentagon_points, 60, 0.4)))
			{
				pent_set1[cnt] = pentagon_points[closest]; 
				closest++;
				cnt++;
			}
			// safety net
			assert (cnt == 5);

			for (cnt=0; cnt < 4; cnt++)
			{
				for (j = cnt + 1; j < 5; )
				{
					if (SquaredDistance (pent_set1[cnt], pent_set1[j]) > 0.4) j++;					
					else 
					{
						if (j != cnt+1)
						{
							temp = pent_set1[cnt+1];
							pent_set1[cnt+1] = pent_set1[j];
							pent_set1[j] = temp;
						}
						break;
					}
				}
			} // end of for()

			index = 0;
			for (cnt = 0; cnt < 5; cnt++)
			{
				for (j = 0; j < seams_parts; j++)
				{
					FindIntermediatePoint (pent_set1[cnt], pent_set1[(cnt+1)%5], polyline[index], ((double)j)/seams_parts);
					index++;
				}
			}

			cnt = 0;
			int positive_start = -1;

			for (index=0; index < total_parts; index++)
			{
				if (GetZCoord (polyline[index]) >= 0) 
				{
					cnt++;
				}
			}

			if (cnt < total_parts && cnt > 0) // not all points are visible
			{
				double z1, z2;

				for (index = 0; index < total_parts; index++)
				{
					if (GetZCoord (polyline[index]) < 0 && GetZCoord (polyline[(index+1)%total_parts]) >= 0)
					{
						positive_start = (index+1) % total_parts;
						z1 = GetZCoord (polyline[index]);
						z2 = GetZCoord (polyline[(index+1)%total_parts]);

						edge[0].Assign (
							polyline[index].x + ((-z1)/(z2-z1))*(polyline[(index+1)%total_parts].x - polyline[index].x),
							polyline[index].y +	((-z1)/(z2-z1))*(polyline[(index+1)%total_parts].y - polyline[index].y),
							polyline[index].z + ((-z1)/(z2-z1))*(polyline[(index+1)%total_parts].z - polyline[index].z));
					}

					if (GetZCoord (polyline[index])<0 && GetZCoord (polyline[(index+total_parts-1)%total_parts])>=0)
					{
						z1 = GetZCoord (polyline[index]);
						z2 = GetZCoord (polyline[(index+total_parts-1)%total_parts]);

						edge[1].Assign (
							polyline[index].x + ((-z1)/(z2-z1))*(polyline[(index+total_parts-1)%total_parts].x - polyline[index].x),
							polyline[index].y +	((-z1)/(z2-z1))*(polyline[(index+total_parts-1)%total_parts].y - polyline[index].y),
							polyline[index].z + ((-z1)/(z2-z1))*(polyline[(index+total_parts-1)%total_parts].z - polyline[index].z));
					}
				}

				index = 0;
				int new_set_cnt = 0;
				index = positive_start;

				double angle = AngleDPoints (edge[0], edge[1]);
				int extra_points = 1 + (int)(angle * 180 / M_PI) / 2;

				DPoint * new_set = new DPoint[total_parts + extra_points + 2];

				new_set[0].Assign (edge[0].x, edge[0].y, edge[0].z);

				while (new_set_cnt < cnt)
				{
					if (GetZCoord (polyline[index])>=0)
					{
						new_set[new_set_cnt + 1].Assign (polyline[index].x,
									polyline[index].y, polyline[index].z);
						new_set_cnt++;
					}

					index++;
					if (index >= total_parts) 
						index = 0;
				}

				new_set[new_set_cnt+1].Assign (edge[1].x, edge[1].y, edge[1].z);
				new_set_cnt++;

				for (index = 1; index < extra_points; index++)
				{
					new_set[new_set_cnt+1].Assign (0,0,0);
					FindIntermediatePoint (edge[1], edge[0], new_set[new_set_cnt+1], ((double)index)/extra_points);
					new_set_cnt++;
				} 

				GeneralPath polygon ;
				for (index = 0; index < new_set_cnt+1; index++)
				{
					if (index == 0)
						polygon.moveTo((float)(centerX + radius * new_set[index].x), 
									 (float)(centerY + radius * new_set[index].y));
					else
						polygon.lineTo((float)(centerX + radius * new_set[index].x), 
									 (float)(centerY + radius * new_set[index].y));
				}

				polygon.closePath();

				polygon.fill(gc, gg);

				delete[] new_set;

			}
			else if (cnt == total_parts)
			{
				GeneralPath polygon ;
				for (index = 0; index < total_parts; index++)
				{
					if (index == 0)
						polygon.moveTo((float)(centerX + radius * polyline[index].x), 
									 (float)(centerY + radius * polyline[index].y));
					else
						polygon.lineTo((float)(centerX + radius * polyline[index].x), 
									 (float)(centerY + radius * polyline[index].y));
				}

				polygon.closePath();

				polygon.fill (gc, gg);

			}	
		}

		delete[] pent_set1;
		delete[] polyline;
		delete[] edge;
	}

	void drawSeams (const DPoint * pentagon_points, GC & gc)
	{
		int j, cnt, index;
		float x, y;
		DPoint * line = new DPoint[seams_parts + 1];
		DPoint * series1 = new DPoint[seams_parts + 1];
		DPoint * series2 = new DPoint[seams_parts + 1];
		DPoint * common = new DPoint [(seams_parts + 1) * 2];
		DPoint cross1, cross2, vector;

		for (int i = 0; i < 60; i++)
		{
			int closest = i;
			while (0 <= (closest = findClosest (pentagon_points[i], closest+1, pentagon_points, 60, 0.4)))
			{
				line[0].Assign (pentagon_points[i].x, pentagon_points[i].y,
										pentagon_points[i].z);

				line[seams_parts].Assign (pentagon_points[closest].x, pentagon_points[closest].y,
										pentagon_points[closest].z);

				// vector is cross product which is perpendicular to line.
				// we use it to create polygon which is 1/50 thick comapred to length of the line.

				vector.Assign (pentagon_points[i].y * pentagon_points[closest].z -
								     pentagon_points[closest].y * pentagon_points[i].z,

									 pentagon_points[i].z * pentagon_points[closest].x -
									 pentagon_points[closest].z * pentagon_points[i].x,

									 pentagon_points[i].x * pentagon_points[closest].y -
									 pentagon_points[closest].x * pentagon_points[i].y);

				vector.Normalize();
				vector.x = vector.x * SEAM_THICKNESS;
				vector.y = vector.y * SEAM_THICKNESS;
				vector.z = vector.z * SEAM_THICKNESS;

				for (j = 1; j <= seams_parts - 1; j++)
				{
					line[j].Assign (0, 0, 0);
					FindIntermediatePoint (pentagon_points[i], pentagon_points[closest], line[j], ((double)j)/seams_parts);
				}

				for (j = 0; j <= seams_parts; j++)
				{
					cross1.x = line[j].x + vector.x;
					cross1.y = line[j].y + vector.y;
					cross1.z = line[j].z + vector.z;

					cross2.x = line[j].x - vector.x;
					cross2.y = line[j].y - vector.y;
					cross2.z = line[j].z - vector.z;

					cross1.Normalize ();
					cross2.Normalize ();

					series1[j].Copy (cross1);
					series2[j].Copy (cross2);
				}

				cnt = 0;

				for (j = 0; j <= seams_parts; j++)
				{
					if (series1[j].z >= 0)
					{
						common [cnt] = series1[j];
						cnt++;
					}
				}

				for (j = seams_parts; j >= 0; j--)
				{
					if (series2[j].z >= 0)
					{
						common[cnt] = series2[j];
						cnt++;
					}
				}

				GeneralPath polygon ;

				index = 0;

				for (j = 0; j < cnt; j++)
				{
					if (common[j].z > 0)
					{
						x = (float)(centerX + radius * common[j].x);
						y =	(float)(centerY + radius * common[j].y);

						if (index == 0)						
							polygon.moveTo (x,y);							
						else
							polygon.lineTo (x, y);

						index++;
					}
				}

				if (index > 1)
				{
					polygon.fill (gc, gg);
				}
			}
		}

		delete[] line;
		delete[] series1;
		delete[] series2;
		delete[] common;
	}

	void FindPentagonPoints (DPoint * pentagon_points)
	{
		int closest;
		int cnt = 0;
		DPoint pent1, pent2;
		
		for (int i=0; i < 12; i++)
		{
			closest = i;
			while (0 <= (closest = findClosest (vertex[i], closest+1, vertex, 12, 1.2)))
			{
				FindIntermediatePoint (vertex[i], vertex[closest], pent1, 1/3.0);
				FindIntermediatePoint (vertex[i], vertex[closest], pent2, 2/3.0);

				if (cnt <= 58)
				{
					pentagon_points[cnt].x = pent1.x;
					pentagon_points[cnt].y = pent1.y;
					pentagon_points[cnt].z = pent1.z;

					cnt++;
					pentagon_points[cnt].x = pent2.x;
					pentagon_points[cnt].y = pent2.y;
					pentagon_points[cnt].z = pent2.z;

					cnt++;
				}
			} 
		}
	}

	void DrawIcosahedron (GC & gc)
	{
		DPoint * pentagon_points = new DPoint[60];

		gg->setPenWidth(gc, 1.0); 

		// this is a little redundant, we can calculate points once
		// instead of doing it on every draw call.
		FindPentagonPoints (pentagon_points);		

		gg->setColor (gc, DKGRAY);

		// "find closest" part in this function also can be done only once.
		drawPaintedPentagons (pentagon_points, gc);

		gg->setColor(gc, BLACK );

		gg->setPenWidth (gc, .62f); 

		drawSeams (pentagon_points, gc);

		gg->drawCircle (gc, centerX, centerY, radius) ;

		delete[] pentagon_points;
	}

	public:


	SoccerBall (int width, int height, GeneralGraphics * _gg) 
	{
		gg = _gg;

		stored = new DPoint[12];
		vertex = new DPoint[12];

		initVertices ();

		radius = BALL_RADIUS ;
		seams_parts = (int)(radius/10);		
		delta = 0;
		tilt = 0;
		delta1 = 0, tilt1 = 0;
		captured = false;
		centerX = width / 2;
		centerY = height / 2;
		iWidth = width;
		iHeight = height;
	}

	~SoccerBall()
	{
		delete[] stored;
		delete[] vertex;
	}

	void resize (int w, int h)
	{
		iWidth = w;
		iHeight = h;

		centerX = w / 2;
		centerY = h / 2;
	}

	void drawAll (GC & gc) 
	{	
		gg->setColor (gc, BACKGD);

		gg->fillRectangle (gc, 0, 0, iWidth, iHeight);

		gg->setColor (gc, WHITE );
		gg->fillCircle (gc, centerX, centerY, radius) ;

		DrawIcosahedron (gc);
	}
};
	
	