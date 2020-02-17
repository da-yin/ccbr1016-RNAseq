#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from lxml import html
import requests
import sys

s = str(sys.argv[1]) 
page = requests.get(s)


tree = html.fromstring(page.content)

# Xpath = '//*[@id="padded_content"]/div[6]/div[1]/div/span'
geneinfo = tree.xpath('//*[@id="padded_content"]/div[6]/div[1]/div/span/text()')

#['Gene ID: 67489, updated on\n                        10-Oct-2019']
try:
	gene = geneinfo[0].split(',')[0].split(':')[-1].strip()
except IndexError:
	gene = 'DOES_NOT_EXIST'

print('{}\t{}'.format(s, gene))


#usage: python xpath_forDa.py 'https://www.ncbi.nlm.nih.gov/gene/?term=ENSMUSG00000099871'
