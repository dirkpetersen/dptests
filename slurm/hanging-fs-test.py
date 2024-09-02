#! /usr/bin/env python3

# Simulate a hanging file system for trouble shooting
#
# Requirements :
#
# dnf install fuse3-devel
# apt install libfuse3-dev
# python3 -m pip install pyfuse3
#
# run as root:
# ./hanging-fs-test.py /mnt/scratch & 
#

import os
import sys
import time
import pyfuse3
import trio
import errno
import stat

class HangingFS(pyfuse3.Operations):
    def __init__(self):
        super().__init__()
        self.files = {}

    async def getattr(self, inode, ctx=None):
        attr = pyfuse3.EntryAttributes()
        attr.st_mode = (stat.S_IFDIR | 0o755)  # Use stat.S_IFDIR for directory type
        attr.st_nlink = 2
        attr.st_size = 0
        attr.st_ctime_ns = int(time.time() * 1e9)
        attr.st_mtime_ns = int(time.time() * 1e9)
        attr.st_atime_ns = int(time.time() * 1e9)
        attr.st_gid = os.getgid()
        attr.st_uid = os.getuid()
        attr.st_ino = inode
        return attr

    async def lookup(self, parent_inode, name, ctx=None):
        raise pyfuse3.FUSEError(errno.ENOENT)

    async def opendir(self, inode, ctx):
        return inode

    async def readdir(self, fh, start_id, token):
        pyfuse3.readdir_reply(token, b'.', await self.getattr(fh), fh)
        pyfuse3.readdir_reply(token, b'..', await self.getattr(fh), fh)

    async def mkdir(self, parent_inode, name, mode, ctx):
        print(f"Hanging on mkdir {name.decode()}")
        while True:
            await trio.sleep(1)

    async def create(self, parent_inode, name, mode, flags, ctx):
        print(f"Hanging on create {name.decode()}")
        while True:
            await trio.sleep(1)

async def main(mountpoint):
    hanging_fs = HangingFS()
    options = set(pyfuse3.default_options)
    options.add('fsname=hanging_fs')
    pyfuse3.init(hanging_fs, mountpoint, options)

    try:
        await pyfuse3.main()
    except:
        pyfuse3.close(unmount=True)
        raise

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <mountpoint>")
        sys.exit(1)

    mountpoint = sys.argv[1]

    trio.run(main, mountpoint)

