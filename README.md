# Sampled ensemble neutrality (SEN)
#

Given a multiple sequence alignment with secondary structure, SEN
calculates the neutrality of the sequences in the alignment and
reports the average. Additionally, SEN can be used to test for
mutational robustness by calculating sequence neutrality and the mean
neutrality of 100 inverse folded sequences, which fold into a similar
structure.

## Usage SEN is written in Clojure 1.5 and requires Java 7+ It is
highly recommended that the user gets Leiningen [here]
(https://github.com/technomancy/leiningen) to run SEN from the command
line.

To calculate the neutrality of sequences in a file
(main-subopt-overlap -f <full-name-of-file>)

Additional help and options
(main-subopt-overlap --help)

If you use the jar file:
java -jar sen-1.0.0-SNAPSHOT-standalone.jar -f "full-name-of-file" --outfile "full-name-of-output"

note the arguments are in double quotes


#TODO
Documentation on API not available yet


## License

;; ------------------------------------------------------------------------ ;;
;; Copyright (c) 2014 Trustees of Boston College                            ;;
;;                                                                          ;;
;; Permission is hereby granted, free of charge, to any person obtaining    ;;
;; a copy of this software and associated documentation files (the          ;;
;; "Software"), to deal in the Software without restriction, including      ;;
;; without limitation the rights to use, copy, modify, merge, publish,      ;;
;; distribute, sublicense, and/or sell copies of the Software, and to       ;;
;; permit persons to whom the Software is furnished to do so, subject to    ;;
;; the following conditions:                                                ;;
;;                                                                          ;;
;; The above copyright notice and this permission notice shall be           ;;
;; included in all copies or substantial portions of the Software.          ;;
;;                                                                          ;;
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,          ;;
;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF       ;;
;; MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                    ;;
;; NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE   ;;
;; LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION   ;;
;; OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION    ;;
;; WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.          ;;
;;                                                                          ;;
;; ------------------------------------------------------------------------ ;;
