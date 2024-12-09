#!/bin/bash

TIMESTAMP=$(date)
sed -i "s|<LAST_UPDATE>|Last updated: $TIMESTAMP|g" README.md