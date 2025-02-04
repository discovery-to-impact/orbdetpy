# server.py - RPC server connection management.
# Copyright (C) 2019-2021 University of Texas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from socket import socket
from atexit import register
from os import environ, path
from psutil import Popen, process_iter
from grpc import channel_ready_future, insecure_channel
from orbdetpy import _data_dir, _jar_file

class RemoteServer:
    rpc_server = None
    rpc_channel = None

    @classmethod
    def connect(cls):
        register(RemoteServer.disconnect)
        jar = path.split(_jar_file)[-1]
        rpc_host, rpc_port = "127.0.0.1", "50051"
        for p in process_iter(attrs=["cmdline"]):
            cmdline = p.info.get("cmdline")
            if (cmdline and cmdline[0].endswith("java")):
                index = next((i for i, x in enumerate(cmdline) if x.endswith(jar)), -1)
                if (index > 0):
                    rpc_port = cmdline[index + 1]
                    break
        else:
            sock = socket()
            sock.bind((rpc_host, 0))
            rpc_port = f"{sock.getsockname()[1]}"
            sock.close()

            java_home = environ.get("JAVA_HOME")
            if (java_home and path.isdir(path.join(java_home, "bin"))):
                java_exec = path.join(java_home, "bin", "java")
            else:
                java_exec = "java"
            cls.rpc_server = Popen([java_exec, "-Xmx2G", "-XX:+UseG1GC", "-jar", _jar_file, rpc_port, _data_dir])

        rpc_uri = f"{rpc_host}:{rpc_port}"
        cls.rpc_channel = insecure_channel(rpc_uri, options=[("grpc.max_send_message_length", 2147483647),
                                                             ("grpc.max_receive_message_length", 2147483647)])
        channel_ready_future(cls.rpc_channel).result(timeout=30.0)

    @classmethod
    def channel(cls):
        if (not cls.rpc_channel):
            RemoteServer.connect()
        return(cls.rpc_channel)

    @classmethod
    def disconnect(cls):
        try:
            if (cls.rpc_channel):
                cls.rpc_channel.close()
            if (cls.rpc_server):
                cls.rpc_server.terminate()
        except Exception as exc:
            pass
