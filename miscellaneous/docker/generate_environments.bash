#!/bin/bash -xe

cp -v ../../.pre-commit-config.yaml build/

MESHTOOLS_WHEEL=${1:-"NO_MESHTOOLS_WHEEL"}
PYCGAL_WHEEL=${2:-"NO_PYCGAL_WHEEL"}
REF_SLUG=${3:+":"}${3:-""}

docker build -t registry.gitlab.inria.fr/charms/compass/base-environment${REF_SLUG} --build-arg MESHTOOLS_WHEEL_PATH=$MESHTOOLS_WHEEL --build-arg PYCGAL_WHEEL_PATH=$PYCGAL_WHEEL -f base/Dockerfile .
docker push registry.gitlab.inria.fr/charms/compass/base-environment${REF_SLUG}

# order is important because of environments dependencies
for evt in build run work profiling doc
do
    echo "Building ${evt}-environment${REF_SLUG}"
    docker build -t registry.gitlab.inria.fr/charms/compass/${evt}-environment${REF_SLUG} ${evt}
    docker push registry.gitlab.inria.fr/charms/compass/${evt}-environment${REF_SLUG}
done
