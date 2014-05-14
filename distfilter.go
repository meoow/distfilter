package main

import "flag"
import . "gomiscutils"
import "bufio"
import "strings"
import "os"
import "fmt"
import "errors"

//command line flags
var snp2genedistfile string
var thresholds string

//temp fields info object
type rs2geneinfo struct {
	rsnum    uint32
	geneid   uint32
	genename string
	dist     uint32
	p        uint8
}

func init() {
	flag.StringVar(&snp2genedistfile, "f", "FILE.TXT", "SNP to GENE distence data file.\n\t(RSID\tGENEID\tGENENAME\tDIST\t5/3)")
	flag.StringVar(&thresholds, "n", "500k", "Threshold (Use just number or with \"k\" suffix)")
}

var p5geneMapped = make(map[uint32]struct{}, 500000)
var hask float64 = 1

func main() {
	flag.Parse()
	if flag.NFlag() < 2 {
		flag.PrintDefaults()
		os.Exit(0)
	}
	//parse the threshold number
	thresholds = strings.ToLower(thresholds)
	if strings.HasSuffix(thresholds,"-") {
		Die(errors.New("Threshold is invalid."))
	}
	if strings.HasSuffix(thresholds, "k") {
		hask = 1000
		thresholds = strings.TrimRight(thresholds, "k")
	}
	threshold := uint32(MustParseFloat(thresholds, 32) * hask)

	//read file the first time to gather info
	//about if 5p dist is satisfied by given threshold
	f, e := os.Open(snp2genedistfile)
	Die(e)
	defer f.Close()
	freader := bufio.NewReader(f)
	lines := Readline(freader)
	for line := range lines {
		fields := strings.Split(TrimNewLine(line), "\t")
		info := fieldParser(fields)
		if info.p == 5 && info.dist <= threshold {
			p5geneMapped[info.rsnum] = struct{}{}
		}
	}
	//reset file and read file the second time to output the rs to gene map
	f.Seek(0, 0)
	lines = Readline(freader)
	for line := range lines {
		fields := strings.Split(TrimNewLine(line), "\t")
		info := fieldParser(fields)
		if info.dist <= threshold {
			switch info.p {
			case 3:
				if _, ok := p5geneMapped[info.rsnum]; ok {
					continue
				}
				//os.Stderr.WriteString("3p found\n")
				fallthrough
			case 5:
				fmt.Printf("%d\t%d\t%s\n", info.rsnum, info.geneid, info.genename)
			default:
				panic("p is invalid")
			}
		}
	}
}

//parse input file text line, as it is function just for code compact
func fieldParser(f []string) *rs2geneinfo {
	rsnum := uint32(MustParseUint(f[0], 10, 32))
	geneid := uint32(MustParseUint(f[1], 10, 32))
	genename := f[2]
	dist := uint32(MustParseUint(f[3], 10, 32))
	p := uint8(MustParseUint(f[4], 10, 8))
	return &rs2geneinfo{rsnum, geneid, genename, dist, p}
}
