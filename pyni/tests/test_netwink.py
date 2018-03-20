#!/usr/bin/env python

from unittest.mock import Mock, MagicMock, patch, call
from unittest import TestCase
import os, shutil, tempfile
from pyni import Netwink

class test_Netwink(TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.netwink = Netwink(self.tempdir,name='testnet',fullInit=False)

    def tearDown(self):
        shutil.rmtree(self.tempdir)
        del self.netwink

    def test_init(self):
        self.assertEqual(self.netwink.name,'testnet')
        self.assertEqual(self.netwink.location,self.tempdir)

    def test_check_location(self):
        self.netwink.check_location()

    def test_set_annotation(self):
        self.netwink.check_location()
        self.netwink.set_annotation()

    def test_load_network(self):
        self.netwink.load_network()

    def test_load_network_cosmic(self):
        self.netwink.load_network(cosmicOnly=True)
