""" A file cache to store data files on a per simulation basis. """
import os
import shutil

class FileCache:
    
    def __init__(self, base, name):
        self.base = base
        self.name = name
        self.cachedir = os.path.join(base, name)
        os.makedirs(self.cachedir, exist_ok=True)

    def cached_file(self, filename, src_dir, invalidate=False):
        
        filepath_in_src = os.path.join(src_dir, filename)
        
        # handle .. in paths
        if ".." in filename:
            filename = filename.replace("..", "__subdir__")
            
        filepath_in_cache = os.path.join(self.cachedir, filename)
        if not os.path.exists(filepath_in_cache) or invalidate:
            os.makedirs(os.path.dirname(filepath_in_cache), exist_ok=True)
            shutil.copy2(filepath_in_src, filepath_in_cache)
        
        return filepath_in_cache