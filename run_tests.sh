#!/usr/bin/bash

docker-compose build cpputest >> /dev/null
docker-compose run cpputest make all
docker-compose rm -f cpputest >> /dev/null
