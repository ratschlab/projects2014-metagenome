#!/bin/bash

# displays a message to stderr in red
function echo_err() {
	RED='\033[0;31m'
	NC='\033[0m'
	echo -e "${RED}Error:${NC} $*" 1>&2;
}

# executes the given command
function execute {
    cmd=("$@")
    echo "Executing ${cmd[*]}"

    # execute the command:
    "${cmd[@]}" || exit 1
}

# check that we can find the aspera ssh key
function check_key {
	declare -a KEY_DIRS=("/usr/local/aspera/connect/etc/asperaweb_id_dsa.openssh" "$HOME/.ssh/asperaweb_id_dsa.openssh")

	for KEY_DIR in "${KEY_DIRS[@]}"; do
		if [ -f "$KEY_DIR" ]; then
			echo "Found aspera key in: $KEY_DIR"
			ASPERA_SSH=$KEY_DIR
			return
		fi
	done
	echo_err "Aspera key not found in any of: " "${KEY_DIRS[@]}" "."
	echo_err "Please copy the key from your installation directory to one of the above locations"
	exit 1
}
