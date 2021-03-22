# PolarProtPred

Prediction of Apical/Basolateral localization from amino acid sequence            
Laszlo Dobson                 
RCNS, Institute of Enzymology // EMBL     
2021                      

## INSTALL
1. Install Keras and TensorFlow: https://phoenixnap.com/kb/how-to-install-keras-on-linux
2. Modify program directory/Blast path at the top of Main.py!
3. Download SwissProt fasta file to the PolarProtPred directory
4. Build blast database from SP fasta file
5. Download speindex.txt from UniProt


## USAGE
    python3 Main.py <FASTA> <TOPOLOGY.XML> <DISORDERED.TXT> <LOW_COMPLEXITY.JSON> <OUTPUT> <MODE>

    FASTA
        Input sequence file, must be in fasta format

    TOPOLOGY.XML
        Topology file in CCTOP XML format
        see: http://cctop.enzim.ttk.mta.hu/?_=/documents/direct_interface.html

    DISORDERED.TXT
        Disordered file in IUPred2A txt format
        see: https://iupred2a.elte.hu/help_new

    LOW_COMPLEXITY.JSON
        Low complexity in PlatoLoCo json format
        see: http://platoloco.aei.polsl.pl/#!/api

    OUTPUT
        Output will be written into the user specified file in the following format:

        Indentifier Loc Score

        The output will also appear on standard output.

    MODE
        binary - apical vs basolateral
        categorical - apical vs basolateral vs other

## Citation

If you find PolarProtPred useful, please cite:

**PolarProtPred: Predicting apical and basolateral localization of transmembrane proteins using putative short linear motifs and deep learning**

Laszlo Dobson, Andras Zeke, Levente Szekeres, Tamas Lango, Gabor Tusnady

Manuscript in preparation.


## PolarProtPred home page
PolarProtPred is available online at [http://polarprotpred.ttk.hu](https://polarprotpred.ttk.hu).

## Related pages
- [PolarProtDb](http://polarprotdb.ttk.hu)
- [CCTOP](http://cctop.enzim.ttk.mta.hu)
- [PDBTM](http://pdbtm.enzim.hu)

