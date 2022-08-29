# RNA2DPred
## Description
Respiratory for the prediction of RNA 2D structure using integer programming.<br/>

RNA2DPred is a Python-based program to predict RNA secondary structure using either energy-based method (MFE), probability-based method (MEA), pattern matching method or 2 of them together (bi-objective optimization).<br/>

The main file for structure prediction is RNA2DPred.py, but the other file motif_alignment.py can be run independently as an alignment tool.

## Usage of motif_alignment
The usage of program motif_alignment.py is like below:<br/>
``` python3  motif_alignment.py --sequence [sequence filename]/[sequence] --motif [motif filename]/[motif sequence] --identity [0 - 1] --format [0,1] --output [output filename] --verbose ```

Run ``` python3  motif_alignment.py --help ``` for more details.

## Usage of RNA2DPred
The usage of program RNA2DPred.py is like below:<br/>
``` python3  RNA2DPred.py --sequence [sequence filename]/[sequence] --motif [motif filename] --identity [0 - 1] --format [0,1] --output [output filename] --verbose ```

Run ``` python3  RNA2DPred.py --help ``` for more details.



