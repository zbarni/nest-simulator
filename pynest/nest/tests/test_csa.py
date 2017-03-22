# -*- coding: utf-8 -*-
#
# test_csa.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

"""
CSA tests
"""

import unittest
import nest

from . import compatibility

try:
    import csa
    HAVE_CSA = True
except ImportError:
    HAVE_CSA = False

try:
    import numpy
    HAVE_NUMPY = True
except ImportError:
    HAVE_NUMPY = False

nest.sli_run("statusdict/have_libneurosim ::")
HAVE_LIBNEUROSIM = nest.sli_pop()


@nest.check_stack
@unittest.skipIf(not HAVE_CSA, 'Python CSA package is not available')
@unittest.skipIf(
    not HAVE_LIBNEUROSIM,
    'PyNEST was built without the libneurosim library'
)
class CSATestCase(unittest.TestCase):
    """CSA tests"""

    def test_CSA_OneToOne_idrange(self):
        """One-to-one connectivity with id ranges"""

        nest.ResetKernel()

        n = 4  # number of neurons

        sources = nest.Create("iaf_neuron", n)
        targets = nest.Create("iaf_neuron", n)

        cg = csa.cset(csa.oneToOne)

        nest.CGConnect(sources, targets, cg)

        for i in range(n):
            conns = nest.GetStatus(
                nest.GetConnections([sources[i]]), 'target')
            self.assertEqual(len(conns), 1)
            self.assertEqual(conns[0], targets[i])

            conns = nest.GetStatus(
                nest.GetConnections([targets[i]]), 'target')
            self.assertEqual(len(conns), 0)

    def test_CSA_OneToOne_params(self):
        """One-to-one connectivity"""

        nest.ResetKernel()

        n = 4  # number of neurons

        sources = nest.Create('iaf_psc_alpha', n)
        targets = nest.Create('iaf_psc_alpha', n)

        cs = csa.cset(csa.oneToOne, 10000.0, 1.0)

        nest.CGConnect(sources, targets, cs, {"weight": 0, "delay": 1})

        for i in range(n):
            conns = nest.GetStatus(
                nest.GetConnections([sources[i]]), 'target')
            self.assertEqual(len(conns), 1)
            self.assertEqual(conns[0], targets[i])

            conns = nest.GetStatus(
                nest.GetConnections([targets[i]]), 'target')
            self.assertEqual(len(conns), 0)

    @unittest.skipIf(not HAVE_NUMPY, 'NumPy package is not available')
    def test_CSA_cgnext(self):
        """cgnext"""

        nest.ResetKernel()

        w = 10000.0
        d = 1.0
        cs = csa.cset(csa.oneToOne, w, d)

        nest.sli_push(cs)
        nest.sli_run('dup')
        nest.sli_push(numpy.array([0, 1, 2, 3]))
        nest.sli_push(numpy.array([0, 1, 2, 3]))
        nest.sli_run('cgsetmask')
        nest.sli_run('dup')
        nest.sli_run('cgstart')
        for i in range(4):
            nest.sli_run('dup')
            nest.sli_run('cgnext')
            self.assertEqual(nest.sli_pop(), True)
            self.assertEqual(nest.sli_pop(), d)
            self.assertEqual(nest.sli_pop(), w)
            self.assertEqual(nest.sli_pop(), i)
            self.assertEqual(nest.sli_pop(), i)
        nest.sli_run('cgnext')
        self.assertEqual(nest.sli_pop(), False)


def suite():

    suite = unittest.makeSuite(CSATestCase, 'test')
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())


if __name__ == "__main__":
    run()
