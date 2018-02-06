#!/bin/sh
./jemdoc -c mysite.conf *.jemdoc
rsync *.html yc6n13@ssh.soton.ac.uk:~/public_html
rsync  *.jpg yc6n13@ssh.soton.ac.uk:~/public_html
rsync  *.mp3 yc6n13@ssh.soton.ac.uk:~/public_html
rsync  *.pdf yc6n13@ssh.soton.ac.uk:~/public_html
echo "Update and upload successful!"
