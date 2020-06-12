#pragma once

#define USE_CAIRO

#include <X11/Xlib.h>
#include <vector>

#ifdef USE_CAIRO
#include <cairo.h>
#include <cairo-xlib.h>
#endif


class GeneralGraphics 
{
	private:

	Visual * visual;
	Display * dpy;
    Pixmap pixmap;

#ifdef USE_CAIRO
	cairo_surface_t *surface = NULL;
    cairo_t *cr = NULL;
#endif

	public:

	Pixmap & getPixmap () 
	{
		return pixmap;
	}

	GeneralGraphics (Display *_dpy, Visual *vis)
	{
		dpy = _dpy;
		visual = vis;
		pixmap = 0L;
	}

	~GeneralGraphics  ()
	{
#ifdef USE_CAIRO
		cairo_destroy (cr);
    	cairo_surface_destroy (surface);
#endif 		
		XFreePixmap(dpy, pixmap); 		
	}

	void resize (int width, int height, int depth, Window & win);

	void setColor (GC & gc, unsigned long RGB);

	void setPenWidth(GC & gc, float width) ;

	void fillRectangle (GC & gc, int left, int top, int right, int bottom);

	void fillCircle (GC & gc, int cx, int cy, float radius) ;

	void drawCircle (GC & gc, int cx, int cy, float radius) ;

	void fillPolygon (GC & gc, const std::vector<std::pair<double,double> > & path);

};
