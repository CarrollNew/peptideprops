*Build*

docker build -t pprops .

*Run*

docker run -dti --name pprops_instance -v /var/log/pprops:/var/log/supervisor/ -p 9000:9000 pprops

