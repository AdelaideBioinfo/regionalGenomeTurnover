// extract fills from net files, min ref size is 10
// Reuben
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

type fill struct {
	RefChr   string
	RefStart int
	RefEnd   int
	QueChr   string
	QueStart int
	QueEnd   int
	Strand   string
	id       int
	//Type     string
}

type gap struct {
	RefChr   string
	RefStart int
	RefEnd   int
	QueStart int
	QueEnd   int
	Strand   string
}

type doubleStack struct {
	ref int
	que int
}

type net struct {
	Chr string
	Len int
}

func main() {

	var newFill fill
	var newGap gap
	var chrStack []string
	var netElem net
	var fillStartStack []doubleStack
	var fillEndStack []doubleStack
	var strandStack []string
	var fillPrint fill
	var chainStack []int

	fillFile, err := os.Open(os.Args[1])
	if err != nil {
		log.Fatalf("reading file: %v", err)
	}
	defer fillFile.Close()

	scFill := bufio.NewScanner(fillFile)
	scFill.Split(bufio.ScanLines)

	for scFill.Scan() {

		var refLen int
		line := scFill.Text()
		line = strings.TrimSpace(line)
		words := strings.Split(line, " ")
		switch words[0] {
		case "net":

			//flush the stack before starting the next chr
			for i := len(fillStartStack) - 1; i > -1; i-- {

				// send fill info to print
				fillPrint.RefStart = fillStartStack[len(fillStartStack)-1].ref
				fillPrint.RefEnd = fillEndStack[len(fillEndStack)-1].ref
				fillPrint.QueStart = fillStartStack[len(fillStartStack)-1].que
				fillPrint.QueEnd = fillEndStack[len(fillEndStack)-1].que

				fillPrint.RefChr = netElem.Chr
				fillPrint.QueChr = chrStack[len(chrStack)-1]
				fillPrint.Strand = strandStack[len(strandStack)-1]

				fillPrint.id = chainStack[len(chainStack)-1]

				//	if fillPrint.RefEnd-fillPrint.RefStart > 8 {
				fmt.Println(strings.Trim(fmt.Sprintf("%v", fillPrint), "{}"))
				//	}

				// remove from stack
				fillStartStack = fillStartStack[0 : len(fillStartStack)-1]
				fillEndStack = fillEndStack[0 : len(fillEndStack)-1]

				chrStack = chrStack[0 : len(chrStack)-1]
				strandStack = strandStack[0 : len(strandStack)-1]
				chainStack = chainStack[0 : len(chainStack)-1]

			}
			netElem.Chr = words[1]
			netElem.Len, err = strconv.Atoi(words[2])

		case "gap":

			// get information for next gap
			refLen, err = strconv.Atoi(words[2])
			if refLen <= 9 {
				continue
			}
			newGap.RefStart, err = strconv.Atoi(words[1])
			newGap.RefStart = newGap.RefStart + 1
			newGap.RefEnd = newGap.RefStart + refLen - 1

			refLen, err = strconv.Atoi(words[6])
			newGap.QueStart, err = strconv.Atoi(words[5])
			newGap.QueStart = newGap.QueStart
			newGap.QueEnd = newGap.QueStart + refLen

			newGap.Strand = words[4]
			if newGap.Strand == "-" {
				newGap.QueStart, newGap.QueEnd = newGap.QueEnd+1, newGap.QueStart-1
			}
			newGap.QueStart = newGap.QueStart + 1
			newGap.QueEnd = newGap.QueEnd

			//fmt.Println(newGap, "newFILL")

			// if the gap is totaly left of our stack, we should print it
			// after that we find we have the fill thats interupted

			if len(fillStartStack) > 0 && len(fillEndStack) > 0 {
				for i := len(fillStartStack) - 1; i > -1; i-- {
					if fillEndStack[i].ref < newGap.RefStart {

						// send fill info to print
						fillPrint.RefStart = fillStartStack[len(fillStartStack)-1].ref
						fillPrint.RefEnd = fillEndStack[len(fillEndStack)-1].ref
						fillPrint.QueStart = fillStartStack[len(fillStartStack)-1].que
						fillPrint.QueEnd = fillEndStack[len(fillEndStack)-1].que

						fillPrint.RefChr = netElem.Chr
						fillPrint.QueChr = chrStack[len(chrStack)-1]
						fillPrint.Strand = strandStack[len(strandStack)-1]
						fillPrint.id = chainStack[len(chainStack)-1]

						//	if fillPrint.RefEnd-fillPrint.RefStart > 8 {
						fmt.Println(strings.Trim(fmt.Sprintf("%v", fillPrint), "{}"))
						//	}

						// remove from stack
						fillStartStack = fillStartStack[0 : len(fillStartStack)-1]
						fillEndStack = fillEndStack[0 : len(fillEndStack)-1]
						chrStack = chrStack[0 : len(chrStack)-1]
						strandStack = strandStack[0 : len(strandStack)-1]
						chainStack = chainStack[0 : len(chainStack)-1]
					}
				}
			}

			if len(fillEndStack) > 0 && newGap.RefStart < fillEndStack[len(fillEndStack)-1].ref {

				fillPrint.RefStart = fillStartStack[len(fillStartStack)-1].ref
				fillPrint.RefEnd = newGap.RefStart - 1 //to get it in position
				fillPrint.QueStart = fillStartStack[len(fillStartStack)-1].que
				fillPrint.QueEnd = newGap.QueStart - 1 //to get it in position

				fillPrint.RefChr = netElem.Chr
				fillPrint.QueChr = chrStack[len(chrStack)-1]
				fillPrint.Strand = strandStack[len(strandStack)-1]
				fillPrint.id = chainStack[len(chainStack)-1]

				//	if fillPrint.RefEnd-fillPrint.RefStart > 8 {
				fmt.Println(strings.Trim(fmt.Sprintf("%v", fillPrint), "{}"))
				//	}

				fillStartStack = fillStartStack[0 : len(fillStartStack)-1]
				fillStartStack = append(fillStartStack, doubleStack{newGap.RefEnd + 1, newGap.QueEnd + 1}) // to get it in position
			}

			//now we need to think what a gap interruption looks like
			// mape our gap interuptions need to work recursivly
			// this may be the problem

		case "fill":

			// get new fill information
			// feed back on position
			newFill.RefStart, err = strconv.Atoi(words[1])
			newFill.RefStart = newFill.RefStart + 1
			refLen, err = strconv.Atoi(words[2])
			newFill.RefEnd = newFill.RefStart + refLen - 1
			newFill.QueChr = words[3]
			newFill.Strand = words[4]
			newFill.QueStart, err = strconv.Atoi(words[5])
			newFill.QueStart = newFill.QueStart
			refLen, err = strconv.Atoi(words[6])
			newFill.QueEnd = newFill.QueStart + refLen
			newFill.id, err = strconv.Atoi(words[8])

			if newFill.Strand == "-" {
				newFill.QueStart, newFill.QueEnd = newFill.QueEnd, newFill.QueStart-1
			}

			newFill.QueStart = newFill.QueStart + 1
			newFill.QueEnd = newFill.QueEnd

			//fmt.Println(newFill, "newGAP")

			// here we enter fill this means we got to the end of our last fill non interrupted
			// here we need to think about our position a bit more
			// If we can get to a start point that is beyond our end point we need to print
			if len(fillStartStack) > 0 && len(fillEndStack) > 0 {
				for i := len(fillStartStack) - 1; i > -1; i-- {
					if fillEndStack[i].ref < newFill.RefStart {

						// send fill info to print
						fillPrint.RefStart = fillStartStack[len(fillStartStack)-1].ref
						fillPrint.RefEnd = fillEndStack[len(fillEndStack)-1].ref
						fillPrint.QueStart = fillStartStack[len(fillStartStack)-1].que
						fillPrint.QueEnd = fillEndStack[len(fillEndStack)-1].que

						fillPrint.RefChr = netElem.Chr
						fillPrint.QueChr = chrStack[len(chrStack)-1]
						fillPrint.Strand = strandStack[len(strandStack)-1]
						fillPrint.id = chainStack[len(chainStack)-1]

						//				if fillPrint.RefEnd-fillPrint.RefStart > 8 {
						fmt.Println(strings.Trim(fmt.Sprintf("%v", fillPrint), "{}"))
						//				}

						// remove from stack
						fillStartStack = fillStartStack[0 : len(fillStartStack)-1]
						fillEndStack = fillEndStack[0 : len(fillEndStack)-1]
						chrStack = chrStack[0 : len(chrStack)-1]
						strandStack = strandStack[0 : len(strandStack)-1]
						chainStack = chainStack[0 : len(chainStack)-1]
					}
				}
			}

			// add to stack
			fillStartStack = append(fillStartStack, doubleStack{newFill.RefStart, newFill.QueStart})
			fillEndStack = append(fillEndStack, doubleStack{newFill.RefEnd, newFill.QueEnd})
			// we could add query start and end to these stacks too maybe?
			// the start can be ref and query, if they both update the same way

			chrStack = append(chrStack, newFill.QueChr)
			strandStack = append(strandStack, newFill.Strand)
			chainStack = append(chainStack, newFill.id)

		}
		if err != nil {
			log.Fatalf("value assignment error", err)
		}
		//fmt.Println(fillStartStack, fillEndStack, "fillStacks")

	}

	for i := len(fillStartStack) - 1; i > -1; i-- {

		// send fill info to print
		fillPrint.RefStart = fillStartStack[len(fillStartStack)-1].ref
		fillPrint.RefEnd = fillEndStack[len(fillEndStack)-1].ref
		fillPrint.QueStart = fillStartStack[len(fillStartStack)-1].que
		fillPrint.QueEnd = fillEndStack[len(fillEndStack)-1].que

		fillPrint.RefChr = netElem.Chr
		fillPrint.QueChr = chrStack[len(chrStack)-1]
		fillPrint.Strand = strandStack[len(strandStack)-1]
		fillPrint.id = chainStack[len(chainStack)-1]

		//	if fillPrint.RefEnd-fillPrint.RefStart > 8 {
		fmt.Println(strings.Trim(fmt.Sprintf("%v", fillPrint), "{}"))
		//	}

		// remove from stack
		fillStartStack = fillStartStack[0 : len(fillStartStack)-1]
		fillEndStack = fillEndStack[0 : len(fillEndStack)-1]
		chrStack = chrStack[0 : len(chrStack)-1]
		strandStack = strandStack[0 : len(strandStack)-1]
		chainStack = chainStack[0 : len(chainStack)-1]

	}

}
