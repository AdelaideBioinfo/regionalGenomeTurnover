// extract gaps from net files, min ref size is 10

// need to see what is happening with the stacking

package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

type gap struct {
	RefChr   string
	RefStart int
	RefEnd   int
	QueChr   string
	QueStart int
	QueEnd   int
	Strand   string
	//Type     string
}

type fill struct {
	RefChr   string
	RefStart int
	RefEnd   int
}

type chainID struct {
	id  int
	end int
}

type net struct {
	Chr string
	Len int
}

func main() {

	var newGap gap
	var newFill fill
	var chrStack []string
	var netElem net
	var gapStartStack []int
	var gapEndStack []int
	var quePosStack []int
	var queEndStack []int
	var strandStack []string
	var gapPrint gap
	var chainInfo chainID
	var chainStack []chainID

	gapFile, err := os.Open(os.Args[1])
	if err != nil {
		log.Fatalf("reading file: %v", err)
	}
	defer gapFile.Close()

	scGap := bufio.NewScanner(gapFile)
	scGap.Split(bufio.ScanLines)

	for scGap.Scan() {

		var refLen int
		line := scGap.Text()
		line = strings.TrimSpace(line)
		words := strings.Split(line, " ")
		switch words[0] {
		case "net":

			// flush the stack
			for i := len(gapStartStack) - 1; i > -1; i-- {

				// send gap info to print
				gapPrint.RefStart = gapStartStack[len(gapStartStack)-1]
				gapPrint.RefEnd = gapEndStack[len(gapEndStack)-1]
				gapPrint.RefChr = netElem.Chr
				gapPrint.QueChr = chrStack[len(chrStack)-1]
				gapPrint.Strand = strandStack[len(strandStack)-1]
				gapPrint.QueStart = quePosStack[len(quePosStack)-1]
				gapPrint.QueEnd = queEndStack[len(queEndStack)-1]

				if gapPrint.RefEnd-gapPrint.RefStart > 8 {
					fmt.Println(strings.Trim(fmt.Sprintf("%v", gapPrint), "{}"), chainStack[len(chainStack)-1].id)
				}
				if gapStartStack[len(gapStartStack)-1] > chainStack[len(chainStack)-1].end {
					chainStack = chainStack[0 : len(chainStack)-1]
				}

				// remove from stack
				gapStartStack = gapStartStack[0 : len(gapStartStack)-1]
				gapEndStack = gapEndStack[0 : len(gapEndStack)-1]
				chrStack = chrStack[0 : len(chrStack)-1]
				quePosStack = quePosStack[0 : len(quePosStack)-1]
				queEndStack = queEndStack[0 : len(queEndStack)-1]
				strandStack = strandStack[0 : len(strandStack)-1]

				// keeping chain stack in sync

			}

			// clean chain stack here
			chainStack = nil

			netElem.Chr = words[1]
			netElem.Len, err = strconv.Atoi(words[2])

		case "fill":

			// if len(gapStartStack) > 0 && gapStartStack[len(gapStartStack)-1] > chainStack[len(chainStack)-1].end {
			// 	chainStack = chainStack[0 : len(chainStack)-1]
			// }

			// get information for next fill

			newFill.RefStart, err = strconv.Atoi(words[1])
			newFill.RefStart = newFill.RefStart + 1
			refLen, err = strconv.Atoi(words[2])
			newFill.RefEnd = newFill.RefStart + refLen - 1

			// if the fill is totally left of our stack, we should print it
			// after that we find we have the gap thats interrupted

			if len(gapStartStack) > 0 && len(gapEndStack) > 0 {
				for i := len(gapStartStack) - 1; i > -1; i-- {
					if gapEndStack[i] < newFill.RefStart {

						// send gap info to print
						gapPrint.RefStart = gapStartStack[len(gapStartStack)-1]
						gapPrint.RefEnd = gapEndStack[len(gapEndStack)-1]
						gapPrint.RefChr = netElem.Chr
						gapPrint.QueChr = chrStack[len(chrStack)-1]
						gapPrint.Strand = strandStack[len(strandStack)-1]
						gapPrint.QueStart = quePosStack[len(quePosStack)-1]
						gapPrint.QueEnd = queEndStack[len(queEndStack)-1]

						if gapPrint.RefStart > chainStack[len(chainStack)-1].end {
							chainStack = chainStack[0 : len(chainStack)-1]
						}

						if gapPrint.RefEnd-gapPrint.RefStart > 8 {
							fmt.Println(strings.Trim(fmt.Sprintf("%v", gapPrint), "{}"), chainStack[len(chainStack)-1].id)
						}

						// remove from stack
						gapStartStack = gapStartStack[0 : len(gapStartStack)-1]
						gapEndStack = gapEndStack[0 : len(gapEndStack)-1]
						chrStack = chrStack[0 : len(chrStack)-1]
						quePosStack = quePosStack[0 : len(quePosStack)-1]
						queEndStack = queEndStack[0 : len(queEndStack)-1]
						strandStack = strandStack[0 : len(strandStack)-1]
					}
					if newFill.RefStart > chainStack[len(chainStack)-1].end {
						chainStack = chainStack[0 : len(chainStack)-1]
					}
				}
			}

			if len(gapEndStack) > 0 && newFill.RefStart < gapEndStack[len(gapEndStack)-1] {

				gapPrint.RefStart = gapStartStack[len(gapStartStack)-1]
				gapPrint.RefEnd = newFill.RefStart - 1
				gapPrint.RefChr = netElem.Chr
				gapPrint.QueChr = chrStack[len(chrStack)-1]
				gapPrint.Strand = strandStack[len(strandStack)-1]
				gapPrint.QueStart = quePosStack[len(quePosStack)-1]
				gapPrint.QueEnd = queEndStack[len(queEndStack)-1]

				if gapPrint.RefEnd-gapPrint.RefStart > 8 {
					fmt.Println(strings.Trim(fmt.Sprintf("%v", gapPrint), "{}"), chainStack[len(chainStack)-1].id)
				}

				gapStartStack = gapStartStack[0 : len(gapStartStack)-1]
				gapStartStack = append(gapStartStack, newFill.RefEnd+1)
			}

			// add all the chain info after considering our fills
			if len(chainStack) > 0 && newFill.RefStart > chainStack[len(chainStack)-1].end {
				chainStack = chainStack[0 : len(chainStack)-1]
			}

			chainInfo.id, err = strconv.Atoi(words[8])
			chainInfo.end = newFill.RefEnd

			chainStack = append(chainStack, chainInfo)

			//now we need to think what a fill interruption looks like
			// mape our fill interuptions need to work recursivly
			// this may be the problem

		case "gap":

			// get new gap information
			// feed back on position
			newGap.RefStart, err = strconv.Atoi(words[1])
			newGap.RefStart = newGap.RefStart + 1
			refLen, err = strconv.Atoi(words[2])
			newGap.RefEnd = newGap.RefStart + refLen - 1
			newGap.QueChr = words[3]
			newGap.Strand = words[4]
			newGap.QueStart, err = strconv.Atoi(words[5])
			newGap.QueStart = newGap.QueStart + 1
			refLen, err = strconv.Atoi(words[6])
			newGap.QueEnd = newGap.QueStart + refLen - 1

			// here we enter gap this means we got to the end of our last gap non interrupted
			// here we need to think about our position a bit more
			// If we can get to a start point that is beyond our end point we need to print
			// chances we have to remove chain information here too
			if len(gapStartStack) > 0 && len(gapEndStack) > 0 {
				for i := len(gapStartStack) - 1; i > -1; i-- {

					if gapEndStack[i] < newGap.RefStart {

						// send gap info to print
						gapPrint.RefStart = gapStartStack[len(gapStartStack)-1]
						gapPrint.RefEnd = gapEndStack[len(gapEndStack)-1]
						gapPrint.RefChr = netElem.Chr
						gapPrint.QueChr = chrStack[len(chrStack)-1]
						gapPrint.Strand = strandStack[len(strandStack)-1]
						gapPrint.QueStart = quePosStack[len(quePosStack)-1]
						gapPrint.QueEnd = queEndStack[len(queEndStack)-1]

						if gapPrint.RefStart > chainStack[len(chainStack)-1].end {
							chainStack = chainStack[0 : len(chainStack)-1]
						}

						if gapPrint.RefEnd-gapPrint.RefStart > 8 {
							fmt.Println(strings.Trim(fmt.Sprintf("%v", gapPrint), "{}"), chainStack[len(chainStack)-1].id)
						}

						// remove from stack
						gapStartStack = gapStartStack[0 : len(gapStartStack)-1]
						gapEndStack = gapEndStack[0 : len(gapEndStack)-1]
						chrStack = chrStack[0 : len(chrStack)-1]
						quePosStack = quePosStack[0 : len(quePosStack)-1]
						queEndStack = queEndStack[0 : len(queEndStack)-1]
						strandStack = strandStack[0 : len(strandStack)-1]
					}

					// if we have gone past our current chain, we remove an element from our stack
					if newGap.RefStart > chainStack[len(chainStack)-1].end {
						chainStack = chainStack[0 : len(chainStack)-1]
					}

				}
			}

			// add to stack
			gapStartStack = append(gapStartStack, newGap.RefStart)
			gapEndStack = append(gapEndStack, newGap.RefEnd)
			chrStack = append(chrStack, newGap.QueChr)
			quePosStack = append(quePosStack, newGap.QueStart)
			queEndStack = append(queEndStack, newGap.QueEnd)
			strandStack = append(strandStack, newGap.Strand)

		}
		if err != nil {
			log.Fatalf("value assignment error", err)
		}

	}

	for i := len(gapStartStack) - 1; i > -1; i-- {

		// send gap info to print
		gapPrint.RefStart = gapStartStack[len(gapStartStack)-1]
		gapPrint.RefEnd = gapEndStack[len(gapEndStack)-1]
		gapPrint.RefChr = netElem.Chr
		gapPrint.QueChr = chrStack[len(chrStack)-1]
		gapPrint.Strand = strandStack[len(strandStack)-1]
		gapPrint.QueStart = quePosStack[len(quePosStack)-1]
		gapPrint.QueEnd = queEndStack[len(queEndStack)-1]

		if gapPrint.RefEnd-gapPrint.RefStart > 8 {
			fmt.Println(strings.Trim(fmt.Sprintf("%v", gapPrint), "{}"), chainStack[len(chainStack)-1].id)
		}

		// remove from stack
		if gapStartStack[len(gapStartStack)-1] > chainStack[len(chainStack)-1].end {
			chainStack = chainStack[0 : len(chainStack)-1]
		}
		gapStartStack = gapStartStack[0 : len(gapStartStack)-1]
		gapEndStack = gapEndStack[0 : len(gapEndStack)-1]
		chrStack = chrStack[0 : len(chrStack)-1]
		quePosStack = quePosStack[0 : len(quePosStack)-1]
		queEndStack = queEndStack[0 : len(queEndStack)-1]
		strandStack = strandStack[0 : len(strandStack)-1]

	}

}
