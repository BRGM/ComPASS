#!/bin/bash

# cf. https://docs.gitlab.com/runner/install/linux-manually.html
apt-get update
apt-get install --yes curl
curl -L --output /usr/local/bin/gitlab-runner https://gitlab-runner-downloads.s3.amazonaws.com/latest/binaries/gitlab-runner-linux-amd64
chmod +x /usr/local/bin/gitlab-runner
curl -sSL https://get.docker.com/ | sh
usermod -aG docker ci
useradd --comment 'GitLab Runner' --create-home gitlab-runner --shell /bin/bash
gitlab-runner install --user=gitlab-runner --working-directory=/home/gitlab-runner
gitlab-runner start
# use --docker-privileged flag or set runner as privileged in /etc/gitlab-runner/config.toml
printf '%s' 'https://gitlab.inria.fr/
XXXX_project_token_XXXX
XXXX_runner_name_XXXX
build,test,run,release
docker
alpine:latest' | gitlab-runner register --docker-privileged
# cf. https://development.robinwinslow.uk/2016/06/23/fix-docker-networking-dns/
echo '{
    "dns": ["172.21.8.87", "193.51.196.130", "193.51.196.131"]
}' > /etc/docker/daemon.json
service docker restart
gitlab-runner restart