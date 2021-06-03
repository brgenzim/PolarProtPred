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

>If you find PolarProtPred useful, please cite:
>
>**PolarProtPred: Predicting apical and basolateral localization of transmembrane proteins using putative short linear motifs and deep learning**
>
>Laszlo Dobson, Andras Zeke, Levente Szekeres, Tamas Lango, Gabor Tusnady
>
>Manuscript in preparation.

## TRAIN

Users can retrain PolarProtPred using the training folder. PolarProtPred ensembles multiple CNNs and fully connected NNs with different 
training and validation sets. The list of proteins for each predictor is included in the TrainValidationSets.csv file (also see original research article
especially Supplementary Table 5 to understand it's logic).

1.   Obtain protein sequences from UniProt. Prepare 2 fasta files for each predictor (training and validation) according to TrainValidationSets.csv.
2.   Predict transmembrane segments using CCTOP. http://cctop.enzim.ttk.mta.hu/?_=/documents/direct_interface.html
3.   Predict disordered regions using IUPred2A. https://iupred2a.elte.hu/help_new
4.   Predict Low complexity regions using PlatoLoCo. http://platoloco.aei.polsl.pl/#!/api
5.   Search for homologous proteins using BLAST on SwissProt - save the result in xml format (-m 7)
6.   Download speindex.txt from UniProt
7.   Generate a label file for each predictor, where each line is space/tab separated, and contains the UniProt AC and the label (1=positive, 0=negative). 
     For ternary prediction 0-1-2 label should be used. 
8.   Supplement these files into subdirectories (BLAST, IDR, LCR, TOP).
9.   Edit train_cnn.py, train_motif.py, train_binary.py, train_ternary.py and modify first lines considering the path of the above generated files. 
     These programs have to be used for multiple times to train all predictors for PolarProtPred
     train_cnn.py - train 5-5-5-5 Apical, Basolateral, Plasma and Other predictor
     train_motif.py - train 5-5-5 Apical, Basolateral and Other predictor
     train_binary.py - train 5 predictor using the output of the previously trained models
     train_ternary.py - train 5 predictor using the output of the previously trained models
(10. We added example test files with 1 sequence to all directories, so it is easier to understand the format of the desired inputs. speindex.txt is not included! Users still have to download
      it from UniProt to try the trainers!)


## PolarProtPred home page
PolarProtPred is available online at [http://polarprotpred.ttk.hu](http://polarprotpred.ttk.hu).

## Related pages
- [PolarProtDb](http://polarprotdb.ttk.hu)
- [CCTOP](http://cctop.enzim.ttk.mta.hu)
- [PDBTM](http://pdbtm.enzim.hu)


