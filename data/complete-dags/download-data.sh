URL="https://zenodo.org/records/17123108/files/complete-dags.tar.xz?download=1"
WORKDIR=$(mktemp -d)

curl $URL -o $WORKDIR/data-supplement.tar.xz
unzip $WORKDIR/data-supplement.zip -d $WORKDIR

mv $WORKDIR/complete-dags/* .
