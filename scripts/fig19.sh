#! /bin/bash

SCRIPT_DIR="$( cd "$(dirname "$0")" ; pwd -P )"
BUILD_DIR=${SCRIPT_DIR}/../build
DATA_DIR=${SCRIPT_DIR}/../data
OUTPUT_DIR=${SCRIPT_DIR}/../data/output

CONTAINER="docker run -it --rm -v ${OUTPUT_DIR}:${OUTPUT_DIR} qnzhou/pyrender"

# Step 1: Download data and compile the code

mkdir -p ${DATA_DIR}
if [ ! -f ${DATA_DIR}/filigree.ply ]; then
	wget https://github.com/geometryprocessing/voroffset/releases/download/1.0/filigree.ply -O ${DATA_DIR}/filigree.ply
fi

mkdir -p ${BUILD_DIR}
if [ ! -f ${BUILD_DIR}/offset3d ]; then
	pushd ${BUILD_DIR}
	cmake ..
	make -j 8
	popd
fi

# Step 2: Run offset algorithm

mkdir -p ${OUTPUT_DIR}

${BUILD_DIR}/offset3d ${DATA_DIR}/filigree.ply -o ${OUTPUT_DIR}/morph_dexelized.obj -f -n 256 -r 3 -p 3 -x noop
${BUILD_DIR}/offset3d ${DATA_DIR}/filigree.ply -o ${OUTPUT_DIR}/morph_dilation.obj -f -n 256 -r 3 -p 3 -x dilation
${BUILD_DIR}/offset3d ${DATA_DIR}/filigree.ply -o ${OUTPUT_DIR}/morph_erosion.obj -f -n 256 -r 3 -p 3 -x erosion
${BUILD_DIR}/offset3d ${DATA_DIR}/filigree.ply -o ${OUTPUT_DIR}/morph_closing.obj -f -n 256 -r 3 -p 3 -x closing
${BUILD_DIR}/offset3d ${DATA_DIR}/filigree.ply -o ${OUTPUT_DIR}/morph_opening.obj -f -n 256 -r 3 -p 3 -x opening

# Step 3: Render pretty pictures

cp ${SCRIPT_DIR}/morph.json ${OUTPUT_DIR}
pushd ${OUTPUT_DIR}
${CONTAINER} bash -c ". /usr/local/mitsuba/setpath.sh; mkdir -p ${OUTPUT_DIR}; cd ${OUTPUT_DIR}; render.py --renderer mitsuba -S `pwd`/morph.json;"
popd
