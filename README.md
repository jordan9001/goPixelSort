# goPixelSort
A pixel sorting implementation in go

## Installation
1. Make sure you have go installed
2. Download "PixelSort.go"
3. Use `go build PixelSort.go` to produce an executable

## Command line Usage
`./goPixelSort -i infile.jpg -o outfile.png -d 750 -r`

+ **-i**  input file, a .jpg or .png (required)
+ **-o**  output file, ending in .jpg or .png (required)
+ **-d**  edge sensitivity, a larger number means more distortion
+ **-r**  reverse sort direction
