# download pyrosetta
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	# pyrosetta python 3.8
	curl 'https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python38.linux/PyRosetta4.Release.python38.linux.release-269.tar.bz2'   -H 'Connection: keep-alive'   -H 'Authorization: Basic bGV2aW50aGFsOnBhcmFkb3g='   -H 'Upgrade-Insecure-Requests: 1'   -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko)
	Chrome/86.0.4240.111 Safari/537.36'   -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9'   -H 'Sec-Fetch-Site: same-origin'   -H 'Sec-Fetch-Mode: navigate'   -H 'Sec-Fetch-Dest: document'   -H 'Referer: https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python38.linux/latest.html'   -H 'Accept-Language: en-GB,en-US;q=0.9,en;q=0.8'   --compressed --output ~/pyrosetta38
        # install vina
	wget http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz -o ~/vina.tar.gz
	tar xfvz ~/vina.tar.gz
	ln -s ~/vina/bin/vina /usr/bin/vina
	ln -s ~/vina/bin/vina_split /usr/bin/vina

elif [[ "$OSTYPE" == "darwin"* ]]; then
        # Mac OSX
	curl 'https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python38.mac/PyRosetta4.Release.python38.mac.release-269.tar.bz2' \
  -H 'Connection: keep-alive' \
  -H 'Authorization: Basic bGV2aW50aGFsOnBhcmFkb3g=' \
  -H 'Upgrade-Insecure-Requests: 1' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/86.0.4240.111 Safari/537.36' \
  -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Sec-Fetch-Mode: navigate' \
  -H 'Sec-Fetch-Dest: document' \
  -H 'Referer: https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python38.mac/latest.html' \
  -H 'Accept-Language: en-GB,en-US;q=0.9,en;q=0.8' \
  -H 'Range: bytes=1523712-1523712' \
  -H 'If-Range: "3dbdec97-5b328ad977000"' \
  --compressed  --output ~/pyrosetta38

cd ~/pyrosetta38/setup
pip install .

pip install biopandas ipython nwalign3

