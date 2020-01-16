#!/bin/bash
# script that gets excecuted when a cloud instance is started

# set up SSD disk
sudo mkfs.ext4 -F /dev/nvme0n1
sudo mkdir -p /mnt/disks/ssd
sudo mount /dev/nvme0n1 /mnt/disks/ssd
sudo chmod a+w /mnt/disks/ssd
echo UUID=`sudo blkid -s UUID -o value /dev/disk/by-id/google-local-nvme-ssd-0` /mnt/disks/ssd ext4 discard,defaults,nofail 0 2 | sudo tee -a /etc/fstab

# get the id of the instance from the metadata server and place it into a file
response=$(curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/instance_id -H "Metadata-Flavor: Google" --write-out  %{http_code} --silent --output /dev/null)
if ((response == 200)); then
  curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/instance_id -H "Metadata-Flavor: Google" > /tmp/instance_id
  instance_id=$(cat /tmp/instance_id)
else
  instance_id='-1'
fi
response=$(curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/server_host -H "Metadata-Flavor: Google" --write-out  %{http_code} --silent --output /dev/null)
if ((response == 200)); then
  curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/server_host -H "Metadata-Flavor: Google" > /tmp/server_host
  server_host=$(cat /tmp/server_host)
else
  server_host='127.0.0.1'
fi

metagraph_path="$HOME/projects2014-metagenome/metagraph"
cd "$metagraph_path"
git checkout dd/cloud  # TODO: remove this at some point
git pull origin dd/cloud

mkdir -p build
cd build
cmake ..
make -j metagraph

cd "$metagraph_path/scripts/cloud"
python3 ./client.py --server_host="${server_host}" --client_id="${instance_id}"