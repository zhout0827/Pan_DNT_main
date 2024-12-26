#!/bin/bash

#Defaultvalueforinputloomfile
DEFAULT_INPUT_LOOM="out.loom"

#Helpfunctiontodisplayusageinformation
usage(){
echo"Usage:$0[OPTIONS]"
echo"OPTIONS:"
echo"-i,--input_loom\tSpecifytheinputloomfilepath"
echo"-h,--help\t\tDisplaythishelpmessage"
exit1
}

#Functiontoprocesscommandlineoptions
process_options(){
whilegetopts"i:h"opt;do
case${opt}in
i)INPUT_LOOM="$OPTARG";;
h)usage;;
\?)echo"Invalidoption:-$OPTARG">&2
exit2;;
:)echo"Option-$OPTARGrequiresanargument.">&2
exit2;;
esac
done
}

#Processtheoptions
process_options"$@"

#Definefilepaths(canusethehg38data)
TFS_PATH="/work/zhout/index/pyscenic/hs_hgnc_tfs.txt"
FEATHER_PATHS="/work/zhout/index/pyscenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
TABLE_PATH="/work/zhout/index/pyscenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

#RunPySCenicGRNanalysiswith4threads
pyscenicgrn\
--num_workers4\
--outputgrn.tsv\
--methodgrnboost2\
"${INPUT_LOOM}""${TFS_PATH}"

#RunPySCenicCistargetanalysiswith14threads
pyscenicctx\
grn.tsv"${FEATHER_PATHS}"\
--annotations_fname"${TABLE_PATH}"\
--expression_mtx_fname"${INPUT_LOOM}"\
--mode"dask_multiprocessing"\
--outputctx.csv\
--num_workers14\
--mask_dropouts

#RunPySCenicAUCellanalysiswith10threads
pyscenicaucell\
"${INPUT_LOOM}"\
ctx.csv\
--outputaucell.loom\
--num_workers10