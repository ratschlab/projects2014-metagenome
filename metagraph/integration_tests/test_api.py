import json
import os
import shlex
import time
from subprocess import Popen
from parameterized import parameterized

import requests

from metagraph.client import Client

from base import TestingBase, METAGRAPH, TEST_DATA_DIR


class TestApi(TestingBase):
    graph_name = 'test_graph'

    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        fasta = TEST_DATA_DIR + '/transcripts_100.fa'

        graph_path = cls.tempdir.name + '/graph.dbg'
        annotation_path = cls.tempdir.name + '/annotation.column.annodbg'

        cls.host = '127.0.0.1'
        cls.port = 3456
        cls._build_graph(cls, fasta, graph_path, 6, 'succinct')

        cls._annotate_graph(cls, fasta, graph_path, annotation_path, 'column')
        print(os.listdir(cls.tempdir.name))

        cls.server_process = cls._start_server(cls, graph_path, annotation_path, cls.port)

        print("waiting for the server to start up on " + str(cls.server_process.pid))
        time.sleep(1)

        cls.graph_client = Client()

        cls.graph_client.add_graph(cls.host, cls.port, cls.graph_name)

    def setUp(self) -> None:
        self.raw_post_request = lambda cmd, payload: requests.post(url=f'http://{self.host}:{self.port}/{cmd}', data=payload)

    @classmethod
    def tearDownClass(cls):
        cls.server_process.kill()

    def _start_server(self, graph, annotation, port):
        construct_command = '{exe} server_query -i {graph} -a {annot} --port {port} --address 127.0.0.1 -p {threads}'.format(
            exe=METAGRAPH,
            graph=graph,
            annot=annotation,
            port=port,
            threads=2
        )

        return Popen(shlex.split(construct_command))

    # do various queries
    def test_api_simple_query(self):
        ret = self.graph_client.search_json("CCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG")

        self.assertIn(self.graph_name, ret.keys())

        res_list, _ = ret[self.graph_name]
        self.assertEqual(len(res_list), 1)

        res_obj = res_list[0]['results']
        self.assertEqual(len(res_obj), 98)

        first_res = res_obj[0]

        self.assertEqual(first_res['sampleCount'], 39)

        self.assertTrue(first_res['sampleName'].startswith('ENST00000456328.2|ENSG00000223972.5|'))
        self.assertTrue('properties' not in first_res.keys()) # doesn't have properties, so don't sent them


    def test_api_raw_incomplete_json(self):
        ret = self.raw_post_request('search', '{"FASTA": ">query\\nAATAAAGGTGTGAGATAACCCCAGCGGTGCCAGGATCCGTGCA", "count_labels": true,')

        self.assertEqual(ret.status_code, 400)
        self.assertIn("Bad json received:", ret.json()['error'])

    def test_api_raw_invalid_params(self):
        # TODO: some parametrized testing, making various cases?
        payload = json.dumps({
                    "num_labels": 'not_a_number',
                    "FASTA": "\n".join([">query",
                                        'AATAAAGGTGTGAGATAACCCCAGCGGTGCCAGGATCCGTGCA',
                                        ]),
                    "count_labels": True,
                    "discovery_fraction": 1 / 100,
                    })

        ret = self.raw_post_request('search', payload)

        self.assertEqual(ret.status_code, 400)
        self.assertIn("Value is not convertible to Int.", ret.json()['error'])

    def test_api_raw_missing_params(self):
        payload = json.dumps({
            "num_labels": 100,
            "count_labels": True,
            "discovery_fraction": 1 / 100
        })

        ret = self.raw_post_request('search', payload)

        self.assertEqual(ret.status_code, 400)
        self.assertIn("No input sequences received from client", ret.json()['error'])

    def test_api_raw_invalid_url(self):
        ret = self.raw_post_request('not_valid', {})
        self.assertEqual(ret.status_code, 404)

    def test_api_raw_no_sequence(self):
        # TODO: some parametrized testing, making various cases?
        payload = json.dumps({
                    "FASTA": "\n".join([">query",
                                        'SEQUENCE_NOT_IN_GRAPH',
                                        ]),
                    "count_labels": True,
                    "discovery_fraction": 1 / 100,
                    "num_labels": 1,
                    })
        ret = self.raw_post_request('search', payload)

        self.assertEqual(ret.status_code, 200)
        self.assertEqual(ret.json(), [])


    def test_api_client_column_labels(self):
        ret = self.graph_client.column_labels()

        self.assertIn(self.graph_name, ret.keys())

        label_list = ret[self.graph_name]

        self.assertGreater(len(label_list), 0)
        self.assertTrue(all(l.startswith('ENST') for l in label_list))


    @parameterized.expand([(1,1), (3,1)])
    def test_api_raw_align_sequence(self, repetitions, foo):
        fasta_str = '\n'.join([ f">query{i}\nTCGATCGA" for i in range(repetitions)])

        payload = json.dumps({"FASTA": fasta_str})

        ret = self.raw_post_request('align', payload)

        self.assertEqual(ret.status_code, 200)

        self.assertEqual(len(ret.json()), repetitions)
        expected = {'score': 12, 'seq_description': 'query0', 'sequence': 'TCGATC'}
        self.assertDictEqual(ret.json()[0], expected)

        self.assertListEqual(
            [ret.json()[i]['seq_description'] for i in range(repetitions)],
            [f"query{i}" for i in range(repetitions)]
        )
