package main

import (
	"flag"
	"fmt"
	"image"
	"image/draw"
	"image/jpeg" // register the JPG format with the image package
	"image/png"  // register the PNG format with the image package
	"os"
	"path/filepath"
)

const (
	VALUE_CONST uint = 15
)

type ArgsIn_t struct {
	help    bool
	file    string
	outfile string
	delta   int
	rev     bool
	angle   int
}

var ArgsIn ArgsIn_t

func init() {
	const (
		default_file  string = ""
		usage_file    string = "(required) Source image (.jpg or .png)"
		usage_outfile string = "(required) Output image (.jpg or .png)"
		default_delta int    = 450
		usage_delta   string = "Sensitivity to edges"
		usage_rev     string = "swap the direction"
		default_angle int    = 0
		usage_angle   string = "(Currently unimplemented)"
	)

	flag.BoolVar(&ArgsIn.help, "h", false, "Prints the usage")

	flag.StringVar(&ArgsIn.file, "image", default_file, usage_file)
	flag.StringVar(&ArgsIn.file, "i", default_file, usage_file+" (shorthand)")
	flag.StringVar(&ArgsIn.outfile, "output", default_file, usage_outfile)
	flag.StringVar(&ArgsIn.outfile, "o", default_file, usage_outfile+" (shorthand)")
	flag.IntVar(&ArgsIn.delta, "delta", default_delta, usage_delta)
	flag.IntVar(&ArgsIn.delta, "d", default_delta, usage_delta+" (shorthand)")
	flag.BoolVar(&ArgsIn.rev, "r", false, usage_rev)
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

	// decode our image
	img, _, err := image.Decode(infile)
	if err != nil {
		fmt.Printf("Unknown file type %s\n", ArgsIn.file)
		flag.PrintDefaults()
		os.Exit(1)
	}

	// put it in a RGBA type
	var rgbaimg *image.RGBA
	if iimg, ok := img.(*image.RGBA); ok {
		rgbaimg = iimg
	} else {
		// it isn't already RGBA, so we have to convert it
		rct := img.Bounds()
		rgbaimg = image.NewRGBA(rct)
		draw.Draw(rgbaimg, rct, img, rct.Min, draw.Src)
	}

	// do rotation
	ArgsIn.angle = (ArgsIn.angle + 360) % 360
	// TODO

	// pixel sort
	pxsrtd, err := PixelSort(rgbaimg, ArgsIn.delta, ArgsIn.rev)
	if err != nil {
		fmt.Printf("Unable to pixel sort due to error : %v\n", err)
		os.Exit(1)
	}

	// rotate back
	// TODO

	// output the image
	outfile, err := os.Create(ArgsIn.outfile)
	if err != nil {
		fmt.Printf("Could not create file %s\n", ArgsIn.outfile)
		os.Exit(1)
	}
	defer outfile.Close()

	if ftype == ".jpg" || ftype == ".jpeg" {
		jpeg.Encode(outfile, pxsrtd, nil)
	} else if ftype == ".png" {
		png.Encode(outfile, pxsrtd)
	}
	os.Exit(0)
}

func PixelSort(img *image.RGBA, delta int, reverse bool) (image.Image, error) {
	// iterate across the image
	bounds := img.Bounds()

	for row := bounds.Min.Y; row < bounds.Max.Y; row++ {
		fmt.Printf("row %d / %d\n", row, bounds.Max.Y-bounds.Min.Y)
		yi := img.Stride * row
		// find areas to sort, within delta

		var back, front int
		back = bounds.Min.X
		for front = back + 1; front < bounds.Max.X; front++ {
			// go until we detect an edge, then sort the area
			i2 := yi + (front * 4)
			i1 := i2 - 4
			p1 := img.Pix[i1 : i1+3]
			p2 := img.Pix[i2 : i2+3]
			// check against delta
			pd := colorDist2(p1, p2)
			if uint32(delta) > pd {
				continue
			}
			// sort and step
			SortArea(img, (back*4)+yi, (front*4)+yi, reverse)
			back = front
		}
		// sort the last bit of the row
		SortArea(img, (back*4)+yi, (front*4)+yi, reverse)
	}

	return img, nil
}

func SortArea(img *image.RGBA, iback, ifrontin int, reverse bool) {
	//fmt.Printf("\t%d %d\n", iback/4, ifrontin/4)
	// bubble pixels
	for ifront := ifrontin; ifront > iback; ifront = ifront - 4 {
		//fmt.Printf("\t\t %d - %d\n", iback, ifront)
		for i := iback; i < ifront-4; i = i + 4 {
			i1 := i
			i2 := i1 + 4
			//fmt.Printf("\t\t%d %d\n", i1, i2)
			p1 := img.Pix[i1 : i1+3]
			p2 := img.Pix[i2 : i2+3]
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
