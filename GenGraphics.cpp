#include <math.h>
#include "GenGraphics.h"

void GeneralGraphics::resize (int width, int height, int depth, Window & win)
{
	if (pixmap)
	{
		XFreePixmap(dpy, pixmap);
	}

	pixmap = XCreatePixmap(dpy, win, width, height, depth);

#ifdef USE_CAIRO	
	if (surface)
	{
		cairo_destroy (cr);
		cr = NULL;
    	cairo_surface_destroy (surface);
		surface = NULL;
	}


	surface = cairo_xlib_surface_create (dpy, pixmap, visual,
					 width, height);

	cr = cairo_create (surface);

	cairo_set_antialias (cr, CAIRO_ANTIALIAS_FAST);
#endif  
}

void GeneralGraphics::setColor (GC & gc, unsigned long RGB)
{
#ifdef USE_CAIRO

	int blue = RGB & 0xff;
	int green = (RGB >> 8) & 0xff;
	int red = (RGB >> 16) & 0xff;

	cairo_set_source_rgb (cr, (float)red/255.0, (float)green / 255.0, (float) blue / 255.0);
#else 		
	XSetForeground(dpy, gc, RGB );
#endif 		

}

void GeneralGraphics::setPenWidth(GC & gc, float width) 
{
#ifdef USE_CAIRO
	cairo_set_line_width(cr, width);
#else 		
	XSetLineAttributes(dpy, gc, width, LineSolid, CapButt, JoinMiter);
#endif 		
}

void GeneralGraphics::fillRectangle (GC & gc, int left, int top, int right, int bottom)
{
#ifdef USE_CAIRO		
	cairo_rectangle (cr, left, top, right - left, bottom - top);
	cairo_fill(cr);
#else 
	XFillRectangle(dpy, pixmap, gc, left, top, right - left, bottom - top);
#endif 		
}

void GeneralGraphics::fillCircle (GC & gc, int cx, int cy, float radius) 
{
#ifdef USE_CAIRO	

	cairo_arc (cr, cx, cy, radius, 0.0, 2 * M_PI);
	cairo_fill(cr);

#else 
	XFillArc (dpy, pixmap, gc, cx - radius, cy - radius, 
				radius * 2, radius * 2, 0, 360 * 64);
#endif 				
}

void GeneralGraphics::drawCircle (GC & gc, int cx, int cy, float radius) 
{
#ifdef USE_CAIRO

	cairo_arc (cr, cx, cy, radius, 0.0, 2 * M_PI);
	cairo_stroke (cr);

#else 		
	XDrawArc (dpy, pixmap, gc, cx - radius, cy - radius, 
				radius * 2, radius * 2, 0, 360 * 64);
#endif 				
}

#ifdef USE_CAIRO
void GeneralGraphics::fillPolygon (GC & gc, const std::vector<std::pair<double,double> > & path)
{
	for (size_t i=0; i < path.size(); i++)
	{
		if (i == 0)
		{
			cairo_move_to(cr, path[i].first, path[0].second);
		}
		else 
		{
			cairo_line_to(cr, path[i].first, path[i].second);
		}
	}

	cairo_fill (cr);
}

#else 
void GeneralGraphics::fillPolygon (GC & gc, const std::vector<std::pair<double,double> > & path)
{
	XPoint *pts = new XPoint [path.size()];

	for (int i=0; i < path.size(); i++)
	{
		pts[i].x = round (path[i].first);
		pts[i].y = round (path[i].second);
	}

	XFillPolygon(dpy, pixmap, gc, pts, path.size(), Nonconvex, CoordModeOrigin);

	delete[] pts;
}
#endif 
