package main

import (
	"flag"
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/jpeg" // register the JPG format with the image package
	"image/png"  // register the PNG format with the image package
	"math"
	"math/cmplx"
	"math/rand"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"

	"github.com/disintegration/imaging"
	"github.com/mjibson/go-dsp/fft"
)

const (
	VALUE_CONST uint = 15
)

const (
	BLOB_SHAPE   int = 0
	BOX_SHAPE    int = 1
	SQUARE_SHAPE int = 2
)

type ArgsIn_t struct {
	help      bool
	file      string
	outfile   string
	delta     int
	rev       bool
	randgroup bool
	angle     float64
	scale     float64
	effect    string
	chanhue   string
	chandeg   int
	chanjoin  bool
	dftsat    float64
}

var ArgsIn ArgsIn_t

func init() {
	const (
		default_file    string  = ""
		usage_file      string  = "(required) Source image (.jpg or .png)"
		usage_outfile   string  = "(required) Output image (.jpg or .png)"
		default_delta   int     = 450
		usage_delta     string  = "Sensitivity to edges for sorting effects"
		usage_rev       string  = "swap the direction" // TODO change this into a angle direction for the grad
		usage_rand      string  = "randomize in sorted group"
		default_scale   float64 = 1.0
		usage_scale     string  = "sample every xth pixel (e.g. 1.5 will downsample to 2/3)"
		default_angle   float64 = 0.0
		usage_angle     string  = "Angle for the applied the effect"
		default_effect  string  = "linesort"
		usage_effect    string  = "Selected effect (linesort|blocksort|floodsort|dft|idft|hdft|hidft)"
		default_chanhue string  = "FF0000"
		default_chandeg int     = 360
		usage_chanhue   string  = "Hue to isolate. Use (r|g|b), sets deg automatically"
		usage_chandeg   string  = "Size of hue to isolate (As degrees, eg 120 for just 1/3 of color space)"
		usage_chanjoin  string  = "Don't join isolated channel after effect"
		usage_dftsat    string  = "Value for log scaling half dft"
		default_dftsat  float64 = 65537.0
	)

	flag.BoolVar(&ArgsIn.help, "h", false, "Prints the usage")

	flag.StringVar(&ArgsIn.file, "image", default_file, usage_file)
	flag.StringVar(&ArgsIn.file, "i", default_file, usage_file+" (shorthand)")
	flag.StringVar(&ArgsIn.outfile, "output", default_file, usage_outfile)
	flag.StringVar(&ArgsIn.outfile, "o", default_file, usage_outfile+" (shorthand)")
	flag.IntVar(&ArgsIn.delta, "delta", default_delta, usage_delta)
	flag.IntVar(&ArgsIn.delta, "d", default_delta, usage_delta+" (shorthand)")
	flag.BoolVar(&ArgsIn.rev, "r", false, usage_rev)
	flag.BoolVar(&ArgsIn.randgroup, "rnd", false, usage_rand)
	flag.Float64Var(&ArgsIn.angle, "ang", default_angle, usage_angle)

	flag.Float64Var(&ArgsIn.scale, "smp", default_scale, usage_scale)

	flag.StringVar(&ArgsIn.effect, "effect", default_effect, usage_effect)
	flag.StringVar(&ArgsIn.effect, "e", default_effect, usage_effect+" (shorthand)")

	flag.StringVar(&ArgsIn.chanhue, "ch", default_chanhue, usage_chanhue)
	flag.IntVar(&ArgsIn.chandeg, "cd", default_chandeg, usage_chandeg)
	flag.BoolVar(&ArgsIn.chanjoin, "cj", false, usage_chanjoin)

	flag.Float64Var(&ArgsIn.dftsat, "dftsat", default_dftsat, usage_dftsat)
	//TODO add a delta mask input option
	//TODO add a dry/wet percentage option
}

func main() {
	// parse the input
	flag.Parse()
	if ArgsIn.help || len(ArgsIn.file) == 0 || len(ArgsIn.outfile) == 0 {
		fmt.Printf("Usage:\n")
		flag.PrintDefaults()
		os.Exit(0)
	}

	infile, err := os.Open(ArgsIn.file)
	if err != nil {
		fmt.Printf("Unknown file %s\n", ArgsIn.file)
		flag.PrintDefaults()
		os.Exit(1)
	}
	defer infile.Close()

	ftype := filepath.Ext(ArgsIn.outfile)
	if ftype != ".jpg" && ftype != ".jpeg" && ftype != ".png" {
		fmt.Printf("Unknown output filetype %s\n", ftype)
		os.Exit(1)
	}

	fmt.Printf("Reading %s\n", ArgsIn.file)
	// decode our image
	img, _, err := image.Decode(infile)
	if err != nil {
		fmt.Printf("Unknown file type %s\n", ArgsIn.file)
		flag.PrintDefaults()
		os.Exit(1)
	}

	// put it in a NRGBA type
	var rgbaimg *image.NRGBA
	if iimg, ok := img.(*image.NRGBA); ok {
		rgbaimg = iimg
	} else {
		// it isn't already RGBA, so we have to convert it
		rct := img.Bounds()
		rgbaimg = image.NewNRGBA(rct)
		draw.Draw(rgbaimg, rct, img, rct.Min, draw.Src)
	}

	// Apply isolation steps (separate channel, downrez, rotate)
	if ArgsIn.scale > 1.0 {
		//TODO scale img down
		obound := rgbaimg.Bounds()
		nw := float64(obound.Dx()) / ArgsIn.scale
		nh := float64(obound.Dy()) / ArgsIn.scale
		nbound := image.Rect(0, 0, int(nw), int(nh))
		newimg := image.NewNRGBA(nbound)
		for y := 0; y < nbound.Max.Y; y++ {
			new_yi := y * newimg.Stride
			old_yi := int(float64(y)*ArgsIn.scale) * rgbaimg.Stride
			for x := 0; x < nbound.Max.X; x++ {
				new_i := new_yi + (x * 4)
				old_i := old_yi + (int(float64(x)*ArgsIn.scale) * 4)

				copy(newimg.Pix[new_i:new_i+4], rgbaimg.Pix[old_i:old_i+4])
			}
		}
		fmt.Printf("Scaled image down to %vx%v\n", int(nw), int(nh))
		rgbaimg = newimg
	}

	// do rotation
	origbounds := rgbaimg.Bounds()
	if ArgsIn.angle != 0 {
		rgbaimg = imaging.Rotate(rgbaimg, ArgsIn.angle, color.Transparent)
	}

	// isolate selected channel
	var otherchans *image.NRGBA = nil

	if ArgsIn.chandeg < 0 || ArgsIn.chandeg > 360 {
		fmt.Println("Please specify a channel spread from 0 to 360")
		os.Exit(1)
	}
	if ArgsIn.chanhue == "r" {
		ArgsIn.chanhue = "FF0000"
		ArgsIn.chandeg = 120
	} else if ArgsIn.chanhue == "g" {
		ArgsIn.chanhue = "00FF00"
		ArgsIn.chandeg = 120
	} else if ArgsIn.chanhue == "b" {
		ArgsIn.chanhue = "0000FF"
		ArgsIn.chandeg = 120
	}
	if ArgsIn.chandeg != 360 {
		c, err := hex2col(ArgsIn.chanhue)
		if err != nil {
			fmt.Printf("Unrecognized hex color %q", ArgsIn.chanhue)
			os.Exit(1)
		}

		if c[3] != 0xff {
			fmt.Printf("Warning, ignoring non-opaque alpha component of channel color")
		}

		if c[0] == 0 && c[1] == 0 && c[2] == 0 {
			fmt.Println("Channel color must not be black")
			os.Exit(1)
		}
		if c[0] == 0xff && c[1] == 0xff && c[2] == 0xff {
			fmt.Println("Channel color must not be white")
			os.Exit(1)
		}

		if ArgsIn.chandeg == 120 {
			if (c[0] == 0xff || c[0] == 0) &&
				(c[1] == 0xff || c[1] == 0) &&
				(c[2] == 0xff || c[2] == 0) {

				// fast mask path
				rgbaimg, otherchans = separateSimpleChannels(rgbaimg, c)
			}

		}

		if otherchans == nil {
			//TODO
			fmt.Println("TODO separate out non-simple RGB channels")
			os.Exit(1)
		}
	}

	if ArgsIn.rev && ArgsIn.randgroup {
		fmt.Printf("Warning, randgroup overrides reverse option")
	}

	// apply selected effect
	fmt.Println("Applying Effect")
	var outimg *image.NRGBA
	if strings.HasPrefix(ArgsIn.effect, "line") || strings.HasPrefix(ArgsIn.effect, "pix") {
		outimg, err = PixelSort(rgbaimg, uint32(ArgsIn.delta), ArgsIn.rev, ArgsIn.randgroup)
		if err != nil {
			fmt.Printf("Unable to pixel sort due to error : %v\n", err)
			os.Exit(1)
		}
	} else if strings.HasPrefix(ArgsIn.effect, "flood") {
		outimg, err = FloodSort(rgbaimg, uint32(ArgsIn.delta), ArgsIn.rev, ArgsIn.randgroup, BLOB_SHAPE)
		if err != nil {
			fmt.Printf("Unable to flood sort due to error : %v\n", err)
			os.Exit(1)
		}
	} else if strings.HasPrefix(ArgsIn.effect, "block") || strings.HasPrefix(ArgsIn.effect, "box") {
		outimg, err = FloodSort(rgbaimg, uint32(ArgsIn.delta), ArgsIn.rev, ArgsIn.randgroup, BOX_SHAPE)
		if err != nil {
			fmt.Printf("Unable to block sort due to error : %v\n", err)
			os.Exit(1)
		}
	} else if strings.HasPrefix(ArgsIn.effect, "square") {
		outimg, err = FloodSort(rgbaimg, uint32(ArgsIn.delta), ArgsIn.rev, ArgsIn.randgroup, SQUARE_SHAPE)
		if err != nil {
			fmt.Printf("Unable to square sort due to error : %v\n", err)
			os.Exit(1)
		}
	} else if strings.HasPrefix(ArgsIn.effect, "hdft") {
		outimg, err = FftHalf(rgbaimg, ArgsIn.dftsat, false)
		if err != nil {
			fmt.Printf("Unable to use fft due to error : %v\n", err)
			os.Exit(1)
		}
	} else if strings.HasPrefix(ArgsIn.effect, "hidft") {
		outimg, err = IFftHalf(rgbaimg, ArgsIn.dftsat, false)
		if err != nil {
			fmt.Printf("Unable to use fft due to error : %v\n", err)
			os.Exit(1)
		}
	} else if strings.HasPrefix(ArgsIn.effect, "dft") {
		outimg, err = FftHalf(rgbaimg, ArgsIn.dftsat, true)
		if err != nil {
			fmt.Printf("Unable to use fft due to error : %v\n", err)
			os.Exit(1)
		}
	} else if strings.HasPrefix(ArgsIn.effect, "idft") {
		outimg, err = IFftHalf(rgbaimg, ArgsIn.dftsat, true)
		if err != nil {
			fmt.Printf("Unable to use fft due to error : %v\n", err)
			os.Exit(1)
		}
	} else {
		fmt.Printf("Unknown effect %v\n", ArgsIn.effect)
		os.Exit(1)
	}

	if otherchans != nil && !ArgsIn.chanjoin {
		// join back with the other channels
		outimg = addImages(outimg, otherchans)
	}

	// rotate back
	if ArgsIn.angle != 0 {
		outimg = imaging.Rotate(rgbaimg, float64(-ArgsIn.angle), color.Transparent)
		// cut back to original size
		newbounds := outimg.Bounds()
		hdx := (newbounds.Dx() - origbounds.Dx()) / 2
		hdy := (newbounds.Dy() - origbounds.Dy()) / 2
		bnd := origbounds.Add(image.Point{X: hdx, Y: hdy})
		outimg = outimg.SubImage(bnd).(*image.NRGBA)
	}

	// output the image
	fmt.Printf("Writing to %s\n", ArgsIn.outfile)
	outfile, err := os.Create(ArgsIn.outfile)
	if err != nil {
		fmt.Printf("Could not create file %s\n", ArgsIn.outfile)
		os.Exit(1)
	}
	defer outfile.Close()

	if ftype == ".jpg" || ftype == ".jpeg" {
		jpeg.Encode(outfile, outimg, nil)
	} else if ftype == ".png" {
		png.Encode(outfile, outimg)
	}
	os.Exit(0)
}

func PixelSort(img *image.NRGBA, delta uint32, reverse bool, randgroup bool) (*image.NRGBA, error) {
	if delta <= 0 {
		return img, nil
	}

	var wg sync.WaitGroup
	// iterate across the image
	bounds := img.Bounds()

	fmt.Printf("Processing %d rows...", bounds.Max.Y-bounds.Min.Y)

	for row := bounds.Min.Y; row < bounds.Max.Y; row++ {
		yi := img.Stride * row
		// find areas to sort, within delta

		wg.Add(1)
		go func() {
			defer wg.Done()
			var back, front int
			back = bounds.Min.X
			inclear := true
			for front = back + 1; front < bounds.Max.X; front++ {

				i2 := yi + (front * 4)
				i1 := i2 - 4

				// first get to a spot past 0 alpha segment
				if img.Pix[i1+3] == 0 {
					back = front
					continue
				}
				inclear = false

				p1 := img.Pix[i1 : i1+4] // front -1
				p2 := img.Pix[i2 : i2+4] // front

				// go until we detect an edge, then sort the area
				dosort := false

				// check if we hit a 0 alpha
				if img.Pix[i2+3] == 0 {
					inclear = true
					dosort = true
				}

				// check against delta
				if !dosort {
					pd := colorDist2(p1, p2)
					if delta < pd {
						dosort = true
					}
				}

				if !dosort {
					continue
				}

				// sort and step
				if !randgroup {
					SortLineArea(img, (back*4)+yi, (front*4)+yi, reverse)
				} else {
					RandomizeLineArea(img, (back*4)+yi, (front*4)+yi)
				}
				back = front
			}
			// sort the last bit of the row as long as it isn't 0 alpha
			if !inclear {
				if !randgroup {
					SortLineArea(img, (back*4)+yi, (front*4)+yi, reverse)
				} else {
					RandomizeLineArea(img, (back*4)+yi, (front*4)+yi)
				}
			}
		}()
	}

	wg.Wait()
	fmt.Printf("\n")
	return img, nil
}

func SortLineArea(img *image.NRGBA, iback, ifrontin int, reverse bool) {
	if iback+4 >= ifrontin {
		return
	}

	// bubble pixels
	for ifront := ifrontin; ifront > iback; ifront = ifront - 4 {
		for i := iback; i < ifront-4; i = i + 4 {
			i1 := i
			i2 := i1 + 4
			p1 := img.Pix[i1 : i1+4]
			p2 := img.Pix[i2 : i2+4]
			sum1 := colorSum(p1)
			sum2 := colorSum(p2)
			// swap the points if brighter on left
			if (reverse && sum1 < sum2) || (!reverse && sum1 > sum2) {
				t := []uint8{0, 0, 0}
				t[0] = p2[0]
				t[1] = p2[1]
				t[2] = p2[2]
				img.Pix[i2] = p1[0]
				img.Pix[i2+1] = p1[1]
				img.Pix[i2+2] = p1[2]
				img.Pix[i1] = t[0]
				img.Pix[i1+1] = t[1]
				img.Pix[i1+2] = t[2]
			}
		}
	}
}

func RandomizeLineArea(img *image.NRGBA, iback, ifront int) {
	rand.Shuffle((ifront-iback)/4, func(i, j int) {
		i1 := iback + (i * 4)
		i2 := iback + (j * 4)
		img.Pix[i1], img.Pix[i2] = img.Pix[i2], img.Pix[i1]
		img.Pix[i1+1], img.Pix[i2+1] = img.Pix[i2+1], img.Pix[i1+1]
		img.Pix[i1+2], img.Pix[i2+2] = img.Pix[i2+2], img.Pix[i1+2]
	})
}

type pxpt struct {
	x int
	y int
}

func FloodSort(img *image.NRGBA, delta uint32, reverse, randgroup bool, shape int) (*image.NRGBA, error) {
	if delta <= 0 {
		return img, nil
	}

	var wg sync.WaitGroup
	// iterate across the image
	bounds := img.Bounds()

	w := bounds.Max.X - bounds.Min.X
	h := bounds.Max.Y - bounds.Min.Y

	if w <= 0 || h <= 0 {
		return img, nil
	}

	var checked []bool = make([]bool, w*h)

	var wave []pxpt = make([]pxpt, 0, w*3)

	// start at min
	wave = append(wave, pxpt{bounds.Min.X, bounds.Min.Y})

	prevcap := cap(wave) // TODO tune this? actually using the biggest is a big mistake though

	debuggroupcount := 0

	//var shapefunc func(startpt pxpt, group *[]pxpt, wave *[]pxpt, img *image.NRGBA, checked []bool, delta uint32)
	var shapefunc func(pxpt, *[]pxpt, *[]pxpt, *image.NRGBA, []bool, uint32, int) = nil

	switch shape {
	case BLOB_SHAPE:
		shapefunc = FloodBlob
	case SQUARE_SHAPE:
		shapefunc = FloodBox
	case BOX_SHAPE:
		shapefunc = FloodBox
	default:
		panic("Unrecognized shape for floodsort")
	}

	// have a processing wave going out, goroutines will sort flooded areas
	for len(wave) > 0 {
		var pt pxpt
		pt, wave = wave[0], wave[1:]

		// if this is already checked, continue
		if checked[(pt.x-bounds.Min.X)+((pt.y-bounds.Min.Y)*w)] {
			continue
		}

		checked[(pt.x-bounds.Min.X)+((pt.y-bounds.Min.Y)*w)] = true

		// if this is 0alpha, then add unchecked neighbors to the wave and continue
		i1 := (pt.y * img.Stride) + (pt.x * 4)
		if img.Pix[i1+3] == 0 {
			for _, npt := range [4]pxpt{{0, 1}, {0, -1}, {1, 0}, {-1, 0}} {
				y := pt.y + npt.y
				x := pt.x + npt.x
				iy := y - bounds.Min.Y
				ix := x - bounds.Min.X

				if ix < 0 || ix >= w || iy < 0 || iy >= h {
					continue
				}

				if checked[ix+(iy*w)] {
					continue
				}

				wave = append(wave, pxpt{x: x, y: y})
			}
			continue
		}

		// flood out a group
		var group []pxpt = make([]pxpt, 0, prevcap)
		debuggroupcount += 1

		shapefunc(pt, &group, &wave, img, checked, delta, shape)

		if len(group) <= 1 {
			continue
		}
		// sort the group in a goroutine
		//fmt.Printf("DEBUG Start %v group for %v (wave at %v) %v/%v\n", debuggroupcount, len(group), len(wave), debugcheckcount, debugfullcheckcount)
		wg.Add(1)
		go func(group []pxpt, img *image.NRGBA, rev bool) {
			defer wg.Done()

			// we should be soul owners of the pixels in the group, so no sync needed beyond waitgroup

			// gather colors
			colors := make([][4]uint8, len(group))
			for i, pt := range group {
				i1 := (pt.y * img.Stride) + (pt.x * 4)
				p1 := img.Pix[i1 : i1+4]
				copy(colors[i][:], p1[:4])
			}

			if !randgroup {
				// sort points by direction
				sort.Slice(group, func(i, j int) bool {
					//DEBUG
					res := group[i].x < group[j].x
					if rev {
						return !res
					}
					return res
				})

				// sort colors
				sort.Slice(colors, func(i, j int) bool {
					sum1 := colorSum(colors[i][:])
					sum2 := colorSum(colors[j][:])
					return sum1 < sum2
				})
			} else {
				// randomize points
				rand.Shuffle(len(group), func(i, j int) {
					group[i], group[j] = group[j], group[i]
				})
			}

			// put sorted colors in the points
			for i := 0; i < len(group); i++ {
				i1 := (group[i].y * img.Stride) + (group[i].x * 4)
				copy(img.Pix[i1:i1+4], colors[i][:])
			}

		}(group, img, reverse)
	}

	fmt.Printf("Waiting on %v sorters\n", debuggroupcount)
	wg.Wait()

	return img, nil
}

func FloodBlob(startpt pxpt, group *[]pxpt, wave *[]pxpt, img *image.NRGBA, checked []bool, delta uint32, _ int) {
	bounds := img.Bounds()
	w := bounds.Dx()
	h := bounds.Dy()

	var groupwave []pxpt = make([]pxpt, 0, cap(*group))
	groupwave = append(groupwave, startpt)

	for len(groupwave) > 0 {
		var gpt pxpt
		gpt, groupwave = groupwave[0], groupwave[1:]

		// add this to the group
		*group = append(*group, gpt)

		// check this is checked
		ix := gpt.x - bounds.Min.X
		iy := gpt.y - bounds.Min.Y
		if !checked[ix+(iy*w)] {
			panic("Px in groupwave that was not already checked")
		}

		i1 := (gpt.y * img.Stride) + (gpt.x * 4)
		p1 := img.Pix[i1 : i1+4]

		// if neighbors are unchecked, see if they are within delta
		// add to group, else add to wave
		for dy := -1; dy <= 1; dy++ {
			if dy == 0 || iy+dy < 0 || iy+dy >= h {
				continue
			}

			if checked[ix+((iy+dy)*w)] {
				continue
			}

			i2 := ((gpt.y + dy) * img.Stride) + (gpt.x * 4)
			p2 := img.Pix[i2 : i2+4]

			// if this is 0alpha add to wave
			addtogroup := true

			if p2[3] == 0 {
				addtogroup = false
			}

			if addtogroup {
				// check if within delta
				pd := colorDist2(p1, p2)
				if delta < pd {
					addtogroup = false
				}
			}

			if !addtogroup {
				// add to wave
				*wave = append(*wave, pxpt{x: gpt.x, y: gpt.y + dy})
			} else {
				// add to group and mark checked
				groupwave = append(groupwave, pxpt{x: gpt.x, y: gpt.y + dy})
				checked[ix+((iy+dy)*w)] = true
			}
		}
		for dx := -1; dx <= 1; dx++ {
			if dx == 0 || ix+dx < 0 || ix+dx >= w {
				continue
			}

			if checked[(ix+dx)+(iy*w)] {
				continue
			}

			i2 := (gpt.y * img.Stride) + ((gpt.x + dx) * 4)
			p2 := img.Pix[i2 : i2+4]

			// if this is 0alpha add to wave
			addtogroup := true

			if p2[3] == 0 {
				addtogroup = false
			}

			if addtogroup {
				// check if within delta
				pd := colorDist2(p1, p2)
				if delta < pd {
					addtogroup = false
				}
			}

			if !addtogroup {
				// add to wave
				*wave = append(*wave, pxpt{x: gpt.x + dx, y: gpt.y})
			} else {
				// add to group
				groupwave = append(groupwave, pxpt{x: gpt.x + dx, y: gpt.y})
				checked[(ix+dx)+(iy*w)] = true
			}
		}
	}
}

func FloodBox(startpt pxpt, group *[]pxpt, wave *[]pxpt, img *image.NRGBA, checked []bool, delta uint32, shape int) {
	// from starting point keep expanding out until we hit something checked, 0alpha, or past delta
	*group = append(*group, startpt)

	bounds := img.Bounds()
	w := bounds.Dx()

	up, down, left, right := 0, 0, 0, 0
	blockup, blockdown, blockleft, blockright := false, false, false, false
	for {
		if shape == SQUARE_SHAPE && ((blockup && blockdown) || (blockleft && blockright)) {
			break
		} else if blockup && blockdown && blockright && blockleft {
			break
		}

		// UP
		if !blockup {
			y := startpt.y - (up + 1)
			if y < bounds.Min.Y {
				blockup = true
			} else {
				for x := startpt.x - left; x <= (startpt.x + right); x++ {
					ix := x - bounds.Min.X
					iy := y - bounds.Min.Y

					if checked[ix+(iy*w)] {
						blockup = true
						break
					}

					i1 := (y * img.Stride) + (x * 4)
					p1 := img.Pix[i1 : i1+4]

					// check 0alpha
					if p1[3] == 0 {
						blockup = true
						break
					}

					i2 := ((y + 1) * img.Stride) + (x * 4)
					p2 := img.Pix[i2 : i2+4]

					pd := colorDist2(p1, p2)
					if pd > delta {
						blockup = true
						break
					}
				}
				if !blockup {
					// add those to the group
					for x := startpt.x - left; x <= (startpt.x + right); x++ {
						ix := x - bounds.Min.X
						iy := y - bounds.Min.Y

						checked[ix+(iy*w)] = true

						*group = append(*group, pxpt{x, y})
					}

					up += 1
				}
			}
		}

		// LEFT
		if !blockleft {
			x := startpt.x - (left + 1)
			if x < bounds.Min.X {
				blockleft = true
			} else {
				for y := startpt.y - up; y <= (startpt.y + down); y++ {
					ix := x - bounds.Min.X
					iy := y - bounds.Min.Y

					if checked[ix+(iy*w)] {
						blockleft = true
						break
					}

					i1 := (y * img.Stride) + (x * 4)
					p1 := img.Pix[i1 : i1+4]

					// check 0alpha
					if p1[3] == 0 {
						blockleft = true
						break
					}

					i2 := (y * img.Stride) + ((x + 1) * 4)
					p2 := img.Pix[i2 : i2+4]

					pd := colorDist2(p1, p2)
					if pd > delta {
						blockleft = true
						break
					}
				}
				if !blockleft {
					// add those to the group
					for y := startpt.y - up; y <= (startpt.y + down); y++ {
						ix := x - bounds.Min.X
						iy := y - bounds.Min.Y

						checked[ix+(iy*w)] = true

						*group = append(*group, pxpt{x, y})
					}

					left += 1
				}
			}
		}

		// DOWN
		if !blockdown {
			y := startpt.y + (down + 1)
			if y >= bounds.Max.Y {
				blockdown = true
			} else {
				for x := startpt.x - left; x <= (startpt.x + right); x++ {
					ix := x - bounds.Min.X
					iy := y - bounds.Min.Y

					if checked[ix+(iy*w)] {
						blockdown = true
						break
					}

					i1 := (y * img.Stride) + (x * 4)
					p1 := img.Pix[i1 : i1+4]

					// check 0alpha
					if p1[3] == 0 {
						blockdown = true
						break
					}

					i2 := ((y - 1) * img.Stride) + (x * 4)
					p2 := img.Pix[i2 : i2+4]

					pd := colorDist2(p1, p2)
					if pd > delta {
						blockdown = true
						break
					}
				}
				if !blockdown {
					// add those to the group
					for x := startpt.x - left; x <= (startpt.x + right); x++ {
						ix := x - bounds.Min.X
						iy := y - bounds.Min.Y

						checked[ix+(iy*w)] = true

						*group = append(*group, pxpt{x, y})
					}

					down += 1
				}
			}
		}

		// RIGHT
		if !blockright {
			x := startpt.x + (right + 1)
			if x >= bounds.Max.X {
				blockright = true
			} else {
				for y := startpt.y - up; y <= (startpt.y + down); y++ {
					ix := x - bounds.Min.X
					iy := y - bounds.Min.Y

					if checked[ix+(iy*w)] {
						blockright = true
						break
					}

					i1 := (y * img.Stride) + (x * 4)
					p1 := img.Pix[i1 : i1+4]

					// check 0alpha
					if p1[3] == 0 {
						blockright = true
						break
					}

					i2 := (y * img.Stride) + ((x - 1) * 4)
					p2 := img.Pix[i2 : i2+4]

					pd := colorDist2(p1, p2)
					if pd > delta {
						blockright = true
						break
					}
				}
				if !blockright {
					// add those to the group
					for y := startpt.y - up; y <= (startpt.y + down); y++ {
						ix := x - bounds.Min.X
						iy := y - bounds.Min.Y

						checked[ix+(iy*w)] = true

						*group = append(*group, pxpt{x, y})
					}

					right += 1
				}
			}
		}
	}

	// add all our nonchecked edges to the wave
	//up & down
	for x := startpt.x - (left + 1); x <= (startpt.x + (right + 1)); x++ {
		if x < bounds.Min.X || x >= bounds.Max.X {
			continue
		}

		ix := x - bounds.Min.X

		y := startpt.y - (up + 1)
		if y >= bounds.Min.Y && y < bounds.Max.Y {
			iy := y - bounds.Min.Y

			if !checked[ix+(iy*w)] {
				*wave = append(*wave, pxpt{x, y})
			}
		}

		y = startpt.y + (down + 1)
		if y >= bounds.Min.Y && y < bounds.Max.Y {
			iy := y - bounds.Min.Y

			if !checked[ix+(iy*w)] {
				*wave = append(*wave, pxpt{x, y})
			}
		}
	}

	//left & right
	for y := startpt.y - (up + 1); y <= (startpt.y + (down + 1)); y++ {
		if y < bounds.Min.Y || y >= bounds.Max.Y {
			continue
		}

		iy := y - bounds.Min.Y

		x := startpt.x - (left + 1)
		if x >= bounds.Min.X && x < bounds.Max.X {
			ix := x - bounds.Min.X

			if !checked[ix+(iy*w)] {
				*wave = append(*wave, pxpt{x, y})
			}
		}

		x = startpt.x + (right + 1)
		if x >= bounds.Min.X && x < bounds.Max.X {
			ix := x - bounds.Min.X

			if !checked[ix+(iy*w)] {
				*wave = append(*wave, pxpt{x, y})
			}
		}
	}

}

func FftHalf(img *image.NRGBA, satval float64, twod bool) (*image.NRGBA, error) {
	fmt.Printf("Floatizing\n")
	// just an image in
	r_real, g_real, b_real := floatize_chans(img, false, satval)

	// parallelize
	fft.SetWorkerPoolSize(0)

	var r_cplx [][]complex128
	var g_cplx [][]complex128
	var b_cplx [][]complex128

	fmt.Printf("fft r\n")
	if twod {
		r_cplx = fft.FFT2Real(r_real)
	} else {
		r_cplx = FftOneChan(r_real)
	}

	fmt.Printf("fft g\n")
	if twod {
		g_cplx = fft.FFT2Real(g_real)
	} else {
		g_cplx = FftOneChan(g_real)
	}

	fmt.Printf("fft b\n")
	if twod {
		b_cplx = fft.FFT2Real(b_real)
	} else {
		b_cplx = FftOneChan(b_real)
	}

	fmt.Printf("together\n")
	// get an image with the mag on the left and the phase on the right
	img = from_complex(r_cplx, g_cplx, b_cplx, true, satval, true)

	return img, nil
}

func FftOneChan(real [][]float64) [][]complex128 {
	cplx := make([][]complex128, len(real))
	for i := 0; i < len(real); i++ {
		cplx[i] = fft.FFTReal(real[i])
	}

	return cplx
}

func IFftHalf(img *image.NRGBA, satval float64, twod bool) (*image.NRGBA, error) {
	fmt.Printf("Floatizing\n")
	r_cplx, g_cplx, b_cplx := complexize_chans(img, true, satval)

	// parallelize
	fft.SetWorkerPoolSize(0)

	fmt.Printf("ifft r\n")
	if twod {
		r_cplx = fft.IFFT2(r_cplx)
	} else {
		r_cplx = IFftOneChan(r_cplx)
	}

	fmt.Printf("ifft g\n")
	if twod {
		g_cplx = fft.IFFT2(g_cplx)
	} else {
		g_cplx = IFftOneChan(g_cplx)
	}

	fmt.Printf("ifft b\n")
	if twod {
		b_cplx = fft.IFFT2(b_cplx)
	} else {
		b_cplx = IFftOneChan(b_cplx)
	}

	fmt.Printf("together\n")
	img = from_complex(r_cplx, g_cplx, b_cplx, false, satval, false)

	return img, nil
}

func IFftOneChan(in [][]complex128) [][]complex128 {
	cplx := make([][]complex128, len(in))
	for i := 0; i < len(in); i++ {
		cplx[i] = fft.IFFT(in[i])
	}

	return cplx
}

func FftMess(img *image.NRGBA, drywet float64, magpix_x, magpix_y, phpix_x, phpix_y int) (*image.NRGBA, error) {

	// parallelize
	fft.SetWorkerPoolSize(0)

	fmt.Printf("Floatizing\n")
	r_real, g_real, b_real := floatize_chans(img, false, 0)

	var r_cplx [][]complex128
	var g_cplx [][]complex128
	var b_cplx [][]complex128

	fmt.Printf("fft r\n")
	r_cplx = fft.FFT2Real(r_real)

	fmt.Printf("fft g\n")
	g_cplx = fft.FFT2Real(g_real)

	fmt.Printf("fft b\n")
	b_cplx = fft.FFT2Real(b_real)

	//TODO
	// my favorites have been pixilating the mag, or phase (separately)
	// or edge detect on the mag overlayed at like 50%
	// or swapping phase with a similar but changed image
	// so, params to add here would be pixilation and edge detect modes
	// and a dry/wet percentage for how much to let the original show through
	// and probably an option to use 1D FFT instead of 2D
	// use delta and delta mask to control amount
	//TODO

	fmt.Printf("ifft b\n")
	b_cplx = fft.IFFT2(b_cplx)

	fmt.Printf("ifft g\n")
	g_cplx = fft.IFFT2(g_cplx)

	fmt.Printf("ifft r\n")
	r_cplx = fft.IFFT2(r_cplx)

	img = from_complex(r_cplx, g_cplx, b_cplx, false, 0, false)

	return img, nil
}

func floatize_chans(img *image.NRGBA, unlog bool, satval float64) (r, g, b [][]float64) {
	ymin := img.Rect.Min.Y
	ymax := img.Rect.Max.Y
	xmin := img.Rect.Min.X
	xmax := img.Rect.Max.X

	r = make([][]float64, ymax-ymin)
	g = make([][]float64, ymax-ymin)
	b = make([][]float64, ymax-ymin)
	for i := 0; i < (ymax - ymin); i++ {
		r[i] = make([]float64, xmax-xmin)
		g[i] = make([]float64, xmax-xmin)
		b[i] = make([]float64, xmax-xmin)
	}

	stride := img.Stride

	scale_c := 255.0 / math.Log(1+satval)

	for y := ymin; y < ymax; y++ {
		y_base := y - ymin
		yoff := y_base * stride
		for x := xmin; x < xmax; x++ {
			x_base := (x - xmin)
			off := yoff + (x_base * 4)

			rf := float64(img.Pix[off+0])
			gf := float64(img.Pix[off+1])
			bf := float64(img.Pix[off+2])

			if unlog {
				// apply a inverst logrithmic transform
				// logged = c ln(1 + v)
				// v = e^(logged/c) - 1
				rf = math.Exp(rf/scale_c) - 1
				gf = math.Exp(gf/scale_c) - 1
				bf = math.Exp(bf/scale_c) - 1
			} else {
				// just scale down to 1
				rf = rf / 255.0
				gf = gf / 255.0
				bf = bf / 255.0
			}

			r[y_base][x_base] = rf
			g[y_base][x_base] = gf
			b[y_base][x_base] = bf
		}
	}

	return r, g, b
}

func complexize_chans(img *image.NRGBA, unlog bool, satval float64) (r, g, b [][]complex128) {
	// assumes right half of image is the phase info
	ymin := img.Rect.Min.Y
	ymax := img.Rect.Max.Y
	xmin := img.Rect.Min.X
	dx := img.Rect.Dx() / 2
	dx4 := dx * 4
	xmax := dx + xmin

	r = make([][]complex128, ymax-ymin)
	g = make([][]complex128, ymax-ymin)
	b = make([][]complex128, ymax-ymin)
	for i := 0; i < (ymax - ymin); i++ {
		r[i] = make([]complex128, xmax-xmin)
		g[i] = make([]complex128, xmax-xmin)
		b[i] = make([]complex128, xmax-xmin)
	}

	stride := img.Stride

	scale_c := 255.0 / math.Log(1+satval)

	phase_map := (2 * math.Pi) / 255.0

	for y := ymin; y < ymax; y++ {
		y_base := y - ymin
		yoff := y_base * stride
		for x := xmin; x < xmax; x++ {
			x_base := (x - xmin)
			off := yoff + (x_base * 4)

			rf := float64(img.Pix[off+0])
			gf := float64(img.Pix[off+1])
			bf := float64(img.Pix[off+2])

			if unlog {
				// apply a inverst logrithmic transform
				// logged = c ln(1 + v)
				// v = e^(logged/c) - 1
				rf = math.Exp(rf/scale_c) - 1
				gf = math.Exp(gf/scale_c) - 1
				bf = math.Exp(bf/scale_c) - 1
			} else {
				// just scale down to 1
				rf = rf / 255.0
				gf = gf / 255.0
				bf = bf / 255.0
			}

			rp := float64(img.Pix[off+dx4+0])
			gp := float64(img.Pix[off+dx4+1])
			bp := float64(img.Pix[off+dx4+2])

			// map back to between -pi and pi
			rp = (rp * phase_map) - math.Pi
			gp = (gp * phase_map) - math.Pi
			bp = (bp * phase_map) - math.Pi

			// z = |z|(cos(p) + isin(p))
			r[y_base][x_base] = complex(rf*math.Cos(rp), rf*math.Sin(rp))
			g[y_base][x_base] = complex(gf*math.Cos(gp), gf*math.Sin(gp))
			b[y_base][x_base] = complex(bf*math.Cos(bp), bf*math.Sin(bp))
		}
	}

	return r, g, b

}

func from_complex(r, g, b [][]complex128, logize bool, satval float64, savephase bool) *image.NRGBA {
	ymax := len(r)
	xmax := len(r[0])

	// monitor for a better satval
	maxabs := 0.0

	// get magnitude
	xwidth := xmax
	if savephase {
		xwidth = xmax * 2
	}

	img := image.NewNRGBA(image.Rectangle{Min: image.Point{0, 0}, Max: image.Point{xwidth, ymax}})
	stride := img.Stride

	scale_c := 255.0 / math.Log(1+satval)

	phase_map := 255.0 / (2 * math.Pi)
	xm4 := xmax * 4

	for y := 0; y < ymax; y++ {
		yoff := y * stride
		for x := 0; x < xmax; x++ {
			off := yoff + (x * 4)

			rabs := cmplx.Abs(r[y][x])
			gabs := cmplx.Abs(g[y][x])
			babs := cmplx.Abs(b[y][x])

			if logize {
				if rabs > maxabs {
					maxabs = rabs
				}
				if gabs > maxabs {
					maxabs = gabs
				}
				if babs > maxabs {
					maxabs = babs
				}

				// apply a logrithmic transform to be able to view it
				// we are losing information here
				rabs = scale_c * math.Log(1+rabs)
				gabs = scale_c * math.Log(1+gabs)
				babs = scale_c * math.Log(1+babs)
			} else {
				// just scale 1 to 256
				rabs *= 255.0
				gabs *= 255.0
				babs *= 255.0
			}

			if rabs > 255 {
				rabs = 255
			} else if rabs < 0 {
				rabs = 0
			}
			if gabs > 255 {
				gabs = 255
			} else if gabs < 0 {
				gabs = 0
			}
			if babs > 255 {
				babs = 255
			} else if babs < 0 {
				babs = 0
			}

			img.Pix[off+0] = byte(rabs)
			img.Pix[off+1] = byte(gabs)
			img.Pix[off+2] = byte(babs)
			img.Pix[off+3] = 0xff

			if savephase {
				rp := cmplx.Phase(r[y][x])
				gp := cmplx.Phase(g[y][x])
				bp := cmplx.Phase(b[y][x])

				// map -pi to pi onto 0 to 255
				rp = (rp + math.Pi) * phase_map
				gp = (gp + math.Pi) * phase_map
				bp = (bp + math.Pi) * phase_map

				img.Pix[off+xm4+0] = byte(rp)
				img.Pix[off+xm4+1] = byte(gp)
				img.Pix[off+xm4+2] = byte(bp)
				img.Pix[off+xm4+3] = 0xff
			}
		}
	}

	if logize {
		fmt.Printf("A better dftsat would have been %v\n", maxabs)
	}

	return img
}

func colorSum(p []uint8) int64 {
	val := int64(p[0]) + int64(p[1]) + int64(p[2])
	var color int64
	color = (int64(p[1]) << 16) + (int64(p[0]) << 8) + (int64(p[2]))
	return (val << VALUE_CONST) + color
}

func colorDist2(p1, p2 []uint8) uint32 {
	var rd, gd, bd int
	// get distance
	rd = int(p1[0]) - int(p2[0])
	rd = rd * rd
	gd = int(p1[1]) - int(p2[1])
	gd = gd * gd
	bd = int(p1[2]) - int(p2[2])
	bd = bd * bd

	return uint32(rd + bd + gd)
}

func hex2col(hex string) ([]uint8, error) {
	cs := make([]uint8, 4)
	cs[3] = 0xff

	hex = strings.TrimPrefix(hex, "#")
	hex = strings.TrimPrefix(hex, "0x")

	var err error = nil

	switch len(hex) {
	case 8:
		_, err = fmt.Sscanf(hex[6:], "%02x", &cs[3])
		if err != nil {
			return nil, err
		}
		fallthrough
	case 6:
		_, err = fmt.Sscanf(hex, "%02x%02x%02x", &cs[0], &cs[1], &cs[2])
	case 4:
		_, err = fmt.Sscanf(hex[3:], "%1x", &cs[3])
		if err != nil {
			return nil, err
		}
		cs[3] |= (cs[3] << 4)
		fallthrough
	case 3:
		_, err = fmt.Sscanf(hex, "%1x%1x%1x", &cs[0], &cs[1], &cs[2])
		cs[0] |= (cs[0] << 4)
		cs[1] |= (cs[1] << 4)
		cs[2] |= (cs[2] << 4)
	default:
		return nil, fmt.Errorf("Color %v is not recognized as a valid hex color length", hex)
	}

	return cs, err
}

func separateSimpleChannels(img *image.NRGBA, c []uint8) (*image.NRGBA, *image.NRGBA) {
	r := c[0] == 0xff
	g := c[1] == 0xff
	b := c[2] == 0xff

	bounds := img.Bounds()
	other := image.NewNRGBA(bounds)

	for y := bounds.Min.Y; y < bounds.Max.Y; y++ {
		yi := img.Stride * y
		for x := bounds.Min.X; x < bounds.Max.X; x++ {
			i := yi + (x * 4)

			if !r {
				other.Pix[i] = img.Pix[i]
				img.Pix[i] = 0
			}

			if !g {
				other.Pix[i+1] = img.Pix[i+1]
				img.Pix[i+1] = 0
			}

			if !b {
				other.Pix[i+2] = img.Pix[i+2]
				img.Pix[i+2] = 0
			}

			// share same alpha I guess
			other.Pix[i+3] = img.Pix[i+3]
		}
	}

	return img, other
}

func addImages(a, b *image.NRGBA) *image.NRGBA {
	if a.Rect != b.Rect {
		panic("Tried to add two images of different sizes")
	}

	bounds := a.Rect

	for y := bounds.Min.Y; y < bounds.Max.Y; y++ {
		yi := a.Stride * y
		for x := bounds.Min.X; x < bounds.Max.X; x++ {
			i := yi + (x * 4)

			a.Pix[i] += b.Pix[i]
			a.Pix[i+1] += b.Pix[i+1]
			a.Pix[i+2] += b.Pix[i+2]
			// Don't add alpha
		}
	}

	return a
}
