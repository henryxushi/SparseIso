PATH_BAMTOOLS_LIB=
BAMFILE=/home/henry/SparseIso_V1.0/demo/demo_single_end.bam
SPARSEISO_PATH=/home/henry/SparseIso_V1.0
OUTPUT=/home/henry/SparseIso_V1.0/demo/demo_single_end
NUM_OF_PROCESSES=5

chmod 777 $SPARSEISO_PATH/SparseIso
chmod 777 $SPARSEISO_PATH/tools/samtools
chmod 777 $SPARSEISO_PATH/tools/processsam

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PATH_BAMTOOLS_LIB
$SPARSEISO_PATH/SparseIso -b $BAMFILE -o $OUTPUT -r s -p $NUM_OF_PROCESSES --cempath $SPARSEISO_PATH/tools --sampath $SPARSEISO_PATH/tools

