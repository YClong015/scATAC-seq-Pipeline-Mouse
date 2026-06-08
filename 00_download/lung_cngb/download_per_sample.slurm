#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=DWLD_CNX0739875
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_nefzger
#SBATCH --output=DWLD_%j.out
#SBATCH --error=DWLD_%j.err


FTP_DIR_URL="ftp://ftp2.cngb.org/pub/CNSA/data6/CNP0004399/CNS0809369/CNX0739880/CNR0841136/"
DOWNLOAD_DIR_NAME="CNX0739880"
DATADIR=/QRISdata/Q8448/Mouse_disease_data/Lung


cd $TMPDIR

echo "Downloading"
wget -c -nH -np -r -R "index.html*" --cut-dirs 5 "$FTP_DIR_URL"
echo "Finish"


if [ ! -d "$DOWNLOAD_DIR_NAME" ]; then
    echo "Wrong didn't create $DOWNLOAD_DIR_NAME"
    exit 1
fi


echo "Download finished, stroe in $DATADIR..."
cp -r ${DOWNLOAD_DIR_NAME} ${DATADIR}/
echo "All complete, data in ${DATADIR}/${DOWNLOAD_DIR_NAME} ---"
