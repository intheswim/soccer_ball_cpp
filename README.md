# Soccer Ball

This is C++ port of my Java applet, using `xlib` and `cairo`.
Original description of the method is available here:

[Mathematics of the soccer ball](http://www.hoist-point.com/soccerball.htm)

## Build instructions (for Ubuntu Linux):

(on other Lixus systems use this link: [Cairo install](https://www.cairographics.org/download/) )

Install cairo dev package:

`sudo apt-get install libcairo2-dev`

Then run make:

`make`

Note: In case you dont need smooth graphics `cairo` provides, you can remove
`#USE_CAIRO` define from `GenGraphics.h` and references to `cairo` from the Makefile. It will be using
default X11 graphics, so image will be pixelated.

## Usage

Use mouse (click, move, release) to rotate the ball.

![Screenshot](/soccerball.png)

## License
[MIT](https://choosealicense.com/licenses/mit/)