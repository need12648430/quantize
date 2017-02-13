# Quantize.js
Quantize.js is a tiny (10.1 kB minified) color quantization library. It extracts approximate color palletes from an input image using a variety of different methods.

It's similar to Color Thief, but twice as fast in Chrome with multiple quantization methods supported.

In addition to the following examples, the source code is clean and documented with [Docco](https://jashkenas.github.io/docco/). It can be found [here](https://need12648430.github.io/quantize/docs).

## Examples
### Find image's dominant color

	var image = document.getElementById("inputImage");
	var quantize = new Quantize();
	var color = quantize.getColor(image);
	// color = [R, G, B];

### Extract 8 colors from image

	var image = document.getElementById("inputImage");
	var quantize = new Quantize();
	var colors = quantize.getPalette(image, 8);
	// color = [[R, G, B], [R, G, B], ...];

Note: The array is ordered by the prominance of the color.
