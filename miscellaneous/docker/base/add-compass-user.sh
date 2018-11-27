#!/bin/bash
source /etc/local/sbin/set-compass-ids.sh
useradd --create-home --shell /bin/bash --uid ${COMPASS_UID} ${COMPASS_USER}
mkdir -p /localfs
chown ${COMPASS_USER}:${COMPASS_USER} /localfs