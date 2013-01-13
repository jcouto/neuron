
import tables as tbl
import numpy as np

def makeOutputFilename(prefix='', extension='.out'):
    filename = prefix
    if prefix != '' and prefix[-1] != '_':
        filename = filename + '_'
    now = time.localtime(time.time())
    filename = filename + '%d%02d%02d-%02d%02d%02d' % \
        (now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)
    if extension[0] != '.':
        extension = '.' + extension
    suffix = ''
    k = 0
    while os.path.exists(filename + suffix + extension):
        k = k+1
        suffix = '_%d' % k
    return filename + suffix + extension

class H5File:
    def __init__(self, filename='', mode='w', title='', compressed=True):
        if filename == '':
            filename = makeOutputFilename(extension='.h5')
        # open the file
        self.fid = tbl.openFile(filename, mode, title)
        self.compressed = compressed
        self.groups = {}
        self.tables = {}
        if self.compressed:
            # create a filter for compression
            self.filter = tbl.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)

    def createGroup(self, where, name, title=''):
        if where[-1] != '/':
            where += '/'
        path = where + name
        try:
            self.groups[path] = self.fid.createGroup(where, name, title)
        except:
            pass

    def writeTable(self, groupName, tableName, description, title, elements):
        if groupName in self.groups:
            self.tables[tableName] = self.fid.createTable(self.groups[groupName], tableName, description, title)
            for k,v in elements.iteritems():
                self.tables[tableName].row[k] = v
            self.tables[tableName].row.append()
            self.tables[tableName].flush()
        else:
            print('Group [%s] not in file.' % groupName)
    
    def writeArray(self, groupName, arrayName, atom, data):
        if type(data) == list:
            data = np.array(data)
        if data.size == 0:
            return
        if groupName in self.groups:
            array = self.fid.createCArray(self.groups[groupName], arrayName, atom, data.shape, filters=self.filter)
            if len(data.shape) == 1:
                array[0:] = data
            else:
                for k in range(data.shape[1]):
                    array[0:,k] = data[0:,k]
        else:
            print('Group [%s] not in file.' % groupName)

    def close(self):
        self.fid.close()



