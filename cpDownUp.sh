

TGTFOLDER=$1

mkdir -p $TGTFOLDER/all
mkdir -p $TGTFOLDER/up
mkdir -p $TGTFOLDER/down

find $TGTFOLDER -maxdepth 1 -type f  | xargs mv -t $TGTFOLDER/all

cp $TGTFOLDER/all/*up*.xlsx $TGTFOLDER/up
cp $TGTFOLDER/all/*up*.pdf $TGTFOLDER/up

cp $TGTFOLDER/all/*down*.xlsx $TGTFOLDER/down
cp $TGTFOLDER/all/*down*.pdf $TGTFOLDER/down

cp /mnt/f/dev/data/Rprojs/covidSC/GO_GSE_RA_README.txt $TGTFOLDER