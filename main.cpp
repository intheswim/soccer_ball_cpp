#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "SoccerBall.h"

struct gameStruct
{
    int width, height, depth;

	GeneralGraphics *gg;

	void cleanup ()
	{
		delete gg;
		gg = nullptr;
	}

	void resize (Window & win, SoccerBall & ball)
	{
		gg->resize (width, height, depth, win);
		ball.resize (width, height);
	}

	void draw (Display *dpy, GC & gc, Window & win)
	{
		XCopyArea(dpy, gg->getPixmap(), win, gc, 0, 0, width, height, 0, 0);

		XFlush (dpy);
	}
} ;


int main (int argc, char *argv[]) 
{
	/* Create our display */
	Display *dpy = XOpenDisplay(NULL);

	/* Get the root window */
	Window root;
	
	/* Let's create our own window. */
	int screen = DefaultScreen(dpy);
	root = XCreateSimpleWindow(dpy, RootWindow(dpy, screen), 24, 48, 860,
				640, 1, BlackPixel(dpy, screen), BlackPixel(dpy, screen));

	XMapWindow(dpy, root);

	XStoreName(dpy, root, "Soccer ball : use mouse to rotate");
	
	XSelectInput (dpy, root, ExposureMask | StructureNotifyMask | ButtonPressMask | ButtonReleaseMask | PointerMotionMask );

	/* Get the window attributes */
	XWindowAttributes wa;
	XGetWindowAttributes(dpy, root, &wa);

	struct gameStruct gs; 

    gs.width = wa.width;
    gs.height = wa.height;
    gs.depth = wa.depth;

	gs.gg = new GeneralGraphics (dpy, DefaultVisual(dpy, screen));

	/* And new create our graphics context to draw on */
	GC gc = XCreateGC(dpy, root, 0, NULL);

	SoccerBall ball (gs.width, gs.height, gs.gg);

	gs.resize (root, ball);

	/* this is to terminate nicely:  */
	Atom wmDeleteMessage = XInternAtom(dpy, "WM_DELETE_WINDOW", False);
    XSetWMProtocols(dpy, root, &wmDeleteMessage, 1);  

	bool keepGoing = true;

	while ( keepGoing )
	{
		XEvent event;

		XNextEvent(dpy, &event);

		switch (event.type)
		{
			case ConfigureNotify: 
        	{
          		XConfigureEvent xce = event.xconfigure;

		        /* This event type is generated for a variety of
          		   happenings, so check whether the window has been
          		   resized. */

          		if (xce.width != gs.width || xce.height != gs.height) 
          		{
            		gs.width = xce.width;
            		gs.height = xce.height;

					gs.resize (root, ball);
          
            		continue;
          		}
        	}
			break;

			case ButtonPress:
			{
				int x = event.xbutton.x;
				int y = event.xbutton.y;

				ball.mouseDown(x, y);
			} 
			break;

			case ButtonRelease:
			{
				int x = event.xbutton.x;
				int y = event.xbutton.y;

				 ball.mouseUp (x, y);
			}
			break;

			case MotionNotify:
			{
				int x = event.xbutton.x;
				int y = event.xbutton.y;

				ball.mouseDrag (x, y);

				XClearArea(dpy, root, 0, 0, 1, 1, true) ; // this causes new Expose message
			}
			break;


			case Expose:
			{
				ball.drawAll(gc);

				gs.draw (dpy, gc, root);
			}
			break;

			case ClientMessage:
        	{
            	if ((Atom)event.xclient.data.l[0] == wmDeleteMessage)
            	{
					keepGoing = false;
                	break;
            	}
			}
			break;

		}

	}

	/* Clear the pixmap used for double buffering */

	gs.cleanup();

	XFreeGC (dpy, gc);
  	XDestroyWindow(dpy, root);
  	XCloseDisplay (dpy);

	return EXIT_SUCCESS;
}
