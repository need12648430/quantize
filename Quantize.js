// # Quantize
// ----------------------------------------------------------------------------
// Quantize is a tiny, but fast color quantization library; it's smaller, and
// twice the speed of ColorThief (with more features, to boot) in Chrome.

/*
Copyright (c) 2017 Jahn Johansen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

var Quantize = (function() {
	"use strict";

	// ## Quantize
	// --------------------------------------------------------------------------
	// This is the Quantize class. It takes an `Object` of settings as an argument.
	// There are simple presets you can pass in defined later. If no settings are
	// supplied, it defaults to a 5-bit modified median cut.
	//
	// Possible settings:
	function Quantize (settings) {
		this.settings = settings || {
			// `histogramPrecision` is a 3-element array defining the precision (in bits)
			// of each channel (R, G, and B).
			histogramPrecision: [5, 5, 5],
			// `fractByPopulace` is used with `Quantize.CutType.Leptonica` to determine
			// when to swich methods.
			// e.g. for 0.75 after 75% of the desired colors are acquired, it will
			// switch from a volume-based cut order to a population-based cut order.
			fractByPopulace: 0.75,
			// `cutType` defines how the bounding box is partitioned. Its values are defined below.
			cutType: CutType.ModifiedMedian,
			// `cutOrder` defines the order in which bounding boxes are partitioned. Its
			// values are defined below.
			cutOrder: CutOrder.Leptonica
		};
	}

	// ### Quantize.CutType
	// --------------------------------------------------------------------------
	// 'enum' of bounding box slicing strategies.
	var CutType = Quantize.CutType = {
		// `Quantize.CutType.Median` is the basic median cut strategy. It cuts along
		// the median of the longest axis (R, G, or B) of the histogram.
		Median: 1,
		// `Quantize.CutType.ModifiedMedian` is a median cut with an offset. Its goal
		// is to better isolate 'unique' colors (rather than averages) by fragmenting
		// less proportionately.
		ModifiedMedian: 2,
		// `Quantize.CutType.Volume` simply cuts the longest axis in half.
		Volume: 3
	};

	// ### Quantize.CutOrder
	// --------------------------------------------------------------------------
	// 'enum' of cut order strategies.
	var CutOrder = Quantize.CutOrder = {
		// `Quantize.CutType.JFIF` is the cut order used to compress standard-compliant
		// JPEG images.
		JFIF: 1,
			// `Quantize.CutType.Leptonica` uses both volume and populace sequentially
			// to decide the partition order.
		Leptonica: 2,
			// `Quantize.CutType.Subdivide` splits *all* boxes until the desired number
			// of colors (or greater) is achieved.
		Subdivide: 3
	};

	// Comparitors
	// --------------------------------------------------------------------------
	// Private. Comparitors used by `Quantize.CutOrder` enums to decide which
	// bounding box to partition next.
	var Comparitor = {
		Population: function(a, b) {
			return a.population() - b.population();
		},
		Volume: function(a, b) {
			return a.volume() - b.volume();
		},
		PopulationVolume: function(a, b) {
			return a.population() * a.volume() - b.population() * b.volume();
		},
		InversePopulation: function(a, b) {
			return b.population() - a.population();
		}
	};

	// ### Quantize.Preset
	// --------------------------------------------------------------------------
	// A collection of quantization strategy presets you can pass to the `Quantize`
	// constructor.
	Quantize.Preset = {
		// `Quantize.Preset.ColorThief` is the Color Thief quantization strategy.
		// It's the modified median cut.
		ColorThief: {
			histogramPrecision: [5, 5, 5],
			fractByPopulace: 0.75,
			cutType: CutType.ModifiedMedian,
			cutOrder: CutOrder.Leptonica
		},
		// `Quantize.Preset.JFIF` is the standard quantization strategy used when
		// compressing JPEG images.
		JFIF: {
			histogramPrecision: [5, 5, 5],
			cutType: CutType.Volume,
			cutOrder: CutOrder.JFIF
		},
		// `Quantize.Preset.MedianCut` is the basic median cut strategy.
		MedianCut: {
			histogramPrecision: [5, 5, 5],
			cutType: CutType.Median,
			cutOrder: CutOrder.Subdivide
		},
		// `Quantize.Preset.Fast` is a modified median cut with a less precise
		// histogram.
		Fast: {
			histogramPrecision: [3, 3, 3],
			fractByPopulace: 0.75,
			cutType: CutType.ModifiedMedian,
			cutOrder: CutOrder.Leptonica
		}
	};

	// ## Histogram
	// --------------------------------------------------------------------------
	// This is the Histogram class used by Quantize. It takes the Quantize
	// instance as its argument to pull its settings.
	function Histogram(quantize) {
		this.precision = quantize.settings.histogramPrecision;

		this.shrink = [
			8 - this.precision[0],
			8 - this.precision[1],
			8 - this.precision[2]
		];
		this.grow = [
			1 << this.shrink[0],
			1 << this.shrink[1],
			1 << this.shrink[2]
		];
		this.mask = [
			(1 << this.precision[0]) - 1,
			(1 << this.precision[1]) - 1,
			(1 << this.precision[2]) - 1
		];

		this.bounds = {
			r: {
				min: 0,
				max: this.mask[0]
			},
			g: {
				min: 0,
				max: this.mask[1]
			},
			b: {
				min: 0,
				max: this.mask[2]
			}
		};

		this.data = new Array(this.mask[0] * this.mask[1] * this.mask[2]);
	}

	Histogram.prototype = {
		// ### Histogram.populate (imageData)
		// ------------------------------------------------------------------------
		// Takes an `ImageData` object as its argument, and populates itself with
		// its colors.
		populate: function(imageData) {
			var bounds = this.bounds, r, g, b, index;

			for (var i = 0; i < imageData.data.length; i += 4) {
				if (imageData.data[i + 3] < 128) continue;

				r = imageData.data[i] >> this.shrink[0];
				g = imageData.data[i + 1] >> this.shrink[1];
				b = imageData.data[i + 2] >> this.shrink[2];

				bounds.r.min = r < bounds.r.min ? r : bounds.r.min;
				bounds.g.min = g < bounds.g.min ? g : bounds.g.min;
				bounds.b.min = b < bounds.b.min ? b : bounds.b.min;
				bounds.r.max = r > bounds.r.max ? r : bounds.r.max;
				bounds.g.max = g > bounds.g.max ? g : bounds.g.max;
				bounds.b.max = b > bounds.b.max ? r : bounds.b.max;

				index = this.index(r, g, b);
				this.data[index] = (this.data[index] || 0) + 1;
			}
		},
		// ### Histogram.index (r, g, b)
		// ------------------------------------------------------------------------
		// A convenience function to map an RGB color to an index in the Histogram
		// instance's internal array.
		index: function(r, g, b) {
			return ((r & this.mask[0]) << (this.precision[0] + this.precision[1])) |
				((g & this.mask[1]) << this.precision[1]) |
				(b & this.mask[2]);
		},
		// ### Histogram.average (bounds)
		// ------------------------------------------------------------------------
		// Calculates the average color within a given bounding box `Box` object.
		average: function(bounds) {
			bounds = typeof bounds === "undefined" ? this.bounds : bounds;

			var r, g, b, value;
			var rsum = 0, bsum = 0, gsum = 0, sum = 0;

			for (r = bounds.r.min; r <= bounds.r.max; r++) {
				for (g = bounds.g.min; g <= bounds.g.max; g++) {
					for (b = bounds.b.min; b <= bounds.b.max; b++) {
						value = this.data[this.index(r, g, b)] || 0;

						rsum += (r + 0.5) * this.grow[0] * value;
						gsum += (g + 0.5) * this.grow[1] * value;
						bsum += (b + 0.5) * this.grow[2] * value;

						sum += value;
					}
				}
			}

			return [~~(rsum / sum), ~~(gsum / sum), ~~(bsum / sum)];
		},
		// ### Histogram.population (bounds)
		// ------------------------------------------------------------------------
		// Calculates the population within a given bounding box `Box` object.
		population: function(bounds) {
			bounds = typeof bounds === "undefined" ? this.bounds : bounds;

			var r, g, b;
			var sum = 0;

			for (r = bounds.r.min; r <= bounds.r.max; r++) {
				for (g = bounds.g.min; g <= bounds.g.max; g++) {
					for (b = bounds.b.min; b <= bounds.b.max; b++) {
						sum += this.data[this.index(r, g, b)] || 0;
					}
				}
			}

			return sum;
		},
		// ### Histogram.population (bounds)
		// ------------------------------------------------------------------------
		// Calculates the smallest bounds that can contain colors in bounding box
		// `Box.`
		trim: function(bounds) {
			bounds = typeof bounds === "undefined" ? this.bounds : bounds;

			var r, g, b, value;

			var trimmed = {
				r: {
					min: this.mask[0] + 1,
					max: -1
				},
				g: {
					min: this.mask[1] + 1,
					max: -1
				},
				b: {
					min: this.mask[2] + 1,
					max: -1
				}
			};

			for (r = bounds.r.min; r <= bounds.r.max; r++) {
				for (g = bounds.g.min; g <= bounds.g.max; g++) {
					for (b = bounds.b.min; b <= bounds.b.max; b++) {
						value = this.data[this.index(r, g, b)] || 0;

						if (value) {
							trimmed.r.min = Math.min(trimmed.r.min, r);
							trimmed.r.max = Math.max(trimmed.r.max, r);
							trimmed.g.min = Math.min(trimmed.g.min, g);
							trimmed.g.max = Math.max(trimmed.g.max, g);
							trimmed.b.min = Math.min(trimmed.b.min, b);
							trimmed.b.max = Math.max(trimmed.b.max, b);
						}
					}
				}
			}

			return trimmed;
		}
	}

	// ## Box
	// --------------------------------------------------------------------------
	// This is the bounding Box class used by Quantize. It represents a subset of the
	// Histogram. It takes the Quantize instance, Histogram instance, and a set of
	// boundaries as its argument.
	//
	// Boundaries are just a JavaScript Object of the form:
	//
	// `{
	// 	r: { min: 0, max: 255 },
	// 	g: { min: 0, max: 255 },
	// 	b: { min: 0, max: 255 }
	// }`
	//
	// It's also responsible for naively memoizing Histogram queries.
	function Box(quantize, histogram, bounds) {
		this.quantize = quantize;
		this.histogram = histogram;
		this.bounds = typeof bounds !== "undefined" ? bounds : histogram.bounds;
		this.memo = {
			populace: null,
			color: null
		}
	}

	Box.prototype = {
		// ### Box.clone ()
		// ------------------------------------------------------------------------
		// Clones this box. Objects pass by reference otherwise.
		clone: function() {
			return new Box(this.quantize, this.histogram, {
				r: {
					min: this.bounds.r.min,
					max: this.bounds.r.max
				},
				g: {
					min: this.bounds.g.min,
					max: this.bounds.g.max
				},
				b: {
					min: this.bounds.b.min,
					max: this.bounds.b.max
				}
			});
		},
		// ### Box.range (axis)
		// ------------------------------------------------------------------------
		// Calculates the length of a given `axis`.
		range: function(axis) {
			return this.bounds[axis].max - this.bounds[axis].min;
		},
		// ### Box.longestAxis ()
		// ------------------------------------------------------------------------
		// Finds the longest axis of this box.
		longestAxis: function() {
			var r = this.range('r'), g = this.range('g'), b = this.range('b');

			if (r >= g && r >= b) return 'r';
			if (g >= r && g >= b) return 'g';
			return 'b';
		},
		// ### Box.color (force)
		// ------------------------------------------------------------------------
		// Calculates the average color of all samples in this Box.
		// If `force` is `true`, it will recalculate a memoized result.
		color: function(force) {
			if (!force && this.memo.color)
				return this.memo.color;

			return this.memo.color = this.histogram.average(this.bounds);
		},
		// ### Box.population ()
		// ------------------------------------------------------------------------
		// Counts and returns all samples contained in this Box.
		// If `force` is `true`, it will recalculate a memoized result.
		population: function(force) {
			if (!force && this.memo.populace)
				return this.memo.populace;

			return this.memo.populace = this.histogram.population(this.bounds);
		},
		// ### Box.volume ()
		// ------------------------------------------------------------------------
		// Calculates the volume of this Box.
		volume: function() {
			return this.range('r') * this.range('g') * this.range('b');
		},
		// ### Box.trim ()
		// ------------------------------------------------------------------------
		// Trims the Box so it's as small as it can be while containing its colors.
		trim: function() {
			this.bounds = this.histogram.trim(this.bounds);
			return this;
		},
		// ### Box.split ()
		// ------------------------------------------------------------------------
		// Splits this box by the `Quantize`-specified cut type.
		split: function() {
			if (this.population() == 1)
				return [this];

			var histogram = this.histogram,
				bounds = this.bounds,
				settings = this.quantize.settings;

			var axis = this.longestAxis(),
				splitPoint;

			if (settings.cutType == CutType.Volume) {
				splitPoint = bounds[axis].min + ~~(this.range(axis) / 2) - 1;
			} else if (settings.cutType == CutType.Median ||
				settings.cutType == CutType.ModifiedMedian) {
				var population = this.population(),
					sum = 0;

				var r = bounds.r,
					g = bounds.g,
					b = bounds.b,
					i, j, k;

				switch (axis) {
					case 'r':
						for (i = r.min; !splitPoint && i <= r.max; i++) {
							for (j = b.min; !splitPoint && j <= b.max; j++) {
								for (k = g.min; !splitPoint && k <= g.max; k++) {
									sum += histogram.data[histogram.index(i, k, j)] || 0;

									if (sum >= population / 2)
										splitPoint = i;
								}
							}
						}
						break;
					case 'g':
						for (i = g.min; !splitPoint && i <= g.max; i++) {
							for (j = b.min; !splitPoint && j <= b.max; j++) {
								for (k = r.min; !splitPoint && k <= r.max; k++) {
									sum += histogram.data[histogram.index(k, i, j)] || 0;

									if (sum >= population / 2)
										splitPoint = i;
								}
							}
						}
						break;
					case 'b':
						for (i = b.min; !splitPoint && i <= b.max; i++) {
							for (j = g.min; !splitPoint && j <= g.max; j++) {
								for (k = r.min; !splitPoint && k <= r.max; k++) {
									sum += histogram.data[histogram.index(k, j, i)] || 0;

									if (sum >= population / 2)
										splitPoint = i;
								}
							}
						}
						break;
				}

				if (settings.cutType == CutType.ModifiedMedian) {
					var left, right;

					left = splitPoint - bounds[axis].min,
						right = bounds[axis].max - splitPoint;

					if (left <= right) {
						splitPoint = Math.min(
							bounds[axis].max - 1, ~~(splitPoint + right / 2)
						);
					} else {
						splitPoint = Math.max(
							bounds[axis].min, ~~(splitPoint - 1 - left / 2)
						);
					}
				}
			}

			var boxA = this.clone(),
				boxB = this.clone();

			boxA.bounds[axis].max = splitPoint;
			boxB.bounds[axis].min = splitPoint + 1;

			return [boxA.trim(), boxB.trim()];
		}
	}

	// ### Quantize.getHistogram (image)
	// ------------------------------------------------------------------------
	// Creates and populates a Histogram from the given Image.
	Quantize.prototype.getHistogram = function (image) {
		var histogram = new Histogram(this);

		var can, ctx, imageData;
		can = document.createElement("canvas");
		can.width = image.width;
		can.height = image.height;
		ctx = can.getContext("2d");
		ctx.drawImage(image, 0, 0);
		imageData = ctx.getImageData(0, 0, image.width, image.height);
		histogram.populate(imageData);

		return histogram;
	}

	// ### Quantize.getPalette (image, targetColors, histogram)
	// ------------------------------------------------------------------------
	// Finds a palette of size `targetColors` that approximates the colors of `image`.
	// `histogram` is optional, if you already populated a Histogram with
	// `Quantize.getHistogram` you can reuse it here.
	Quantize.prototype.getPalette = function (image, targetColors, histogram) {
		targetColors = typeof targetColors === "undefined" ? 4 : targetColors;

		histogram = typeof histogram === "undefined" ? this.getHistogram(image) : histogram;
		var boxes = [new Box(this, histogram)], box, split, lastComparitor;

		function insert(array, value, comparitor, i, tmp) {
			array.push(value);
			if (array.length == 1) return;
			for (i = array.length - 1; i > 0; i--) {
				tmp = array[i];
				if (comparitor(tmp, array[i - 1]) < 0) {
					array[i] = array[i - 1];
					array[i - 1] = tmp;
				}
			}
		}

		function splitNextBox(comparitor) {
			if (comparitor != lastComparitor)
				boxes.sort(lastComparitor = comparitor);

			var box, split;
			box = boxes.pop();

			split = box.split();

			while (split.length)
				insert(boxes, split.pop(), comparitor)
		}

		switch (this.settings.cutOrder) {
			case CutOrder.Leptonica:
				while (boxes.length < targetColors * this.settings.fractByPopulace)
					splitNextBox(Comparitor.Population, true);

				while (boxes.length < targetColors)
					splitNextBox(Comparitor.PopulationVolume);
				break;
			case CutOrder.JFIF:
				while (boxes.length < targetColors / 2)
					splitNextBox(Comparitor.Population, true);

				while (boxes.length < targetColors)
					splitNextBox(Comparitor.Volume);
				break;
			case CutOrder.Subdivide:
				while (boxes.length < targetColors) {
					var amount = boxes.length;
					for (var i = 0; i < amount; i++) {
						split = boxes.shift().split();
						while (split.length)
							boxes.push(split.shift());
					}
				}
				break;
		}

		boxes.sort(Comparitor.InversePopulation);

		for (var i = 0; i < boxes.length; i++)
			boxes[i] = boxes[i].color();

		return boxes;
	}

	// ### Quantize.quantizeImage (image, targetColors)
	// ------------------------------------------------------------------------
	// Lazily quantizes an image. It's not fast, as this isn't actually what I
	// wrote the library to do. I just wanted to extract pretty color palettes.
	Quantize.prototype.quantizeImage = function (image, targetColors) {
		function difference (a, b) {
			return Math.abs(Math.sqrt(
				Math.pow(b[0] - a[0], 2) +
				Math.pow(b[1] - a[1], 2) +
				Math.pow(b[2] - a[2], 2)
			));
		}

		var canvas, ctx, imageData, histogram, palette;
		canvas = document.createElement("canvas");
		canvas.width = image.width;
		canvas.height = image.height;

		ctx = canvas.getContext("2d");
		ctx.drawImage(image, 0, 0);
		imageData = ctx.getImageData(0, 0, image.width, image.height);

		histogram = new Histogram(this);
		histogram.populate(imageData);
		palette = this.getPalette(image, targetColors, histogram);

		var i, j, pixel, distance, closestIndex, closestDistance;

		for (i = 0; i < imageData.data.length; i += 4) {
			pixel = [imageData.data[i], imageData.data[i + 1], imageData.data[i + 2]]
			closestIndex = 0;
			closestDistance = 256;

			for (j = 0; j < palette.length; j ++) {
				distance = difference(pixel, palette[j]);
				if (distance < closestDistance) {
					closestDistance = distance;
					closestIndex = j;
				}
			}

			imageData.data[i] = palette[closestIndex][0];
			imageData.data[i + 1] = palette[closestIndex][1];
			imageData.data[i + 2] = palette[closestIndex][2];
		}

		ctx.putImageData(imageData, 0, 0);
		return canvas.toDataURL("image/png");
	}

	// ### Quantize.getColor (image)
	// ------------------------------------------------------------------------
	// Tries to find an images most dominant color, similar to Color Thief.
	Quantize.prototype.getColor = function (image) {
		return this.getPalette(image, 5).shift();
	}

	return Quantize;
})();
