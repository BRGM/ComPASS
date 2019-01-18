source /etc/local/sbin/set-compass-ids.sh
echo ""
echo "This is a docker container made for ComPASS."
if [ -z $COMPASS_UID ]
then
    echo "The ComPASS user UID is not set!"
else
    if [ $UID -eq $COMPASS_UID ]
    then
        echo "Default user is $COMPASS_USER with UID 1000"
        echo "If necessary set the user UID using the --user docker run flag"
    else
        echo "You are running container with UID ${UID} which is different"
        echo "from the ComPASS default user UID which is $COMPASS_UID"
        echo "This is why you won't have user name on this session"
    fi
fi
echo ""
export PS1='ComPASS:\w\$ '
# The following is to pass mpi/petsc flags easily to cmake
export CC=mpicc
