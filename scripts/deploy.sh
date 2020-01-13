#!/usr/bin/env bash
if [[ "$#" == 0 ]]; then
	HOST="localhost"
elif [[ "$#" == 1 ]]; then
	HOST=$1
else
	echo "Wrong syntax! Usage: deploy.sh 'hostname'"
	exit 1
fi

if [[ ! -e "src" || ! -e "setup.py" ]]; then
	echo "Could not find src dir or setup.py. Make sure they exist!"
	exit 1
else
	CODENAME="$(basename $(realpath .))"
fi

UUID="$(uuidgen)"
REMOTE_TMP_DIR="/tmp/$UUID-install-$CODENAME"

echo "Deploying '$CODENAME' on '$HOST'"
if [[ "$HOST" == "localhost" ]]; then
	# copy files
	rsync -r --exclude "src/*.egg-info" src setup.py $REMOTE_TMP_DIR
	# install and clean up
	sh -c "cd $REMOTE_TMP_DIR; python3 setup.py install --user; cd -; rm -rf $REMOTE_TMP_DIR" 1>/dev/null
else
	# copy files
	rsync -r --exclude "src/*.egg-info" src setup.py $HOST:$REMOTE_TMP_DIR
	# install and clean up
	ssh $HOST "cd $REMOTE_TMP_DIR; python3 setup.py install --user; cd -; rm -rf $REMOTE_TMP_DIR" 1>/dev/null
fi
