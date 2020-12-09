#! /usr/bin/env python3

# basic example for reading and writing csv files directly from 
# Google drive into a pandas dataframe. 
# All files are merged into a single data frame 
# and written back to google drive 
# credentials are stored in ~/.gcp 

import pickle, os, pandas, io, time
from googleapiclient import discovery 
from googleapiclient.http import MediaIoBaseDownload
from googleapiclient.http import MediaIoBaseUpload
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request

# If modifying these scopes, delete the file ~/.gcp/token.pickle .
#SCOPES = ['https://www.googleapis.com/auth/drive.metadata.readonly']
SCOPES = ['https://www.googleapis.com/auth/drive']

def main():
    """Shows basic usage of the Drive v3 API.
    Prints the names and ids of the first 10 files the user has access to.
    """
    service = discovery.build('drive', 'v3', credentials=getCredentials())

    # Call the Drive v3 API to get all the csv files in Drive
    results = service.files().list(
        pageSize=100, 
        q="mimeType='text/csv'", 
        spaces='drive',
        fields="nextPageToken, files(id, name, size)").execute()
    items = results.get('files', [])

    bigdf=None
    for item in items:
        print(u'\n{0} ({1}) {2} bytes'.format(item['name'],  item['id'], item['size']))
        if not item['name'].startswith('users'):
            continue  
        t1 = time.time()
        obj = pullObject(service, item['id'])
        df = pandas.read_csv(io.StringIO(obj.decode()))
        sec = time.time() - t1
        MBS =  int(item['size'])/1048576/sec
        print(df)
        print('{0} download speed: {1:.3f} MB/s\n'.format(item['name'],  MBS)) 
        
        # merge data frames 
        if isinstance(bigdf, pandas.DataFrame):
            bigdf =  bigdf.append(df)            
        else:
            bigdf = df
            
    print(bigdf)
    myStream = io.BytesIO()  # io.StringIO() #
    #bigdf.to_csv(myStream.encode())
    #ret=pushObject(service, "moin.csv",  myStream)

def pullObject(service, file_id, progress=True):
    #request = service.files().export(fileId=file_id, mimeType='text/csv')
    request = service.files().get_media(fileId=file_id)
    fh = io.BytesIO()
    downloader = MediaIoBaseDownload(fh, request)
    done = False
    while done is False:
        status, done = downloader.next_chunk()
        if progress:
            print("Download %d%%." % int(status.progress() * 100))
    return fh.getvalue()

def pushObject(service, file_name, content, progress=True):    
    #https://stackoverflow.com/questions/36214138/google-drive-python-api-upload-a-media-object-which-is-not-a-file
    #https://stackoverflow.com/questions/50535181/cant-upload-a-csv-as-google-sheet-with-special-characters-in-the-contents

    #with open('MyFile.jpeg', 'rb') as FID:
    #     fileInMemory = FID.read()

    media = MediaIoBaseUpload(content, mimetype='application/octet-stream', resumable=False) #io.BytesIO(content) mimetype='text/csv' 
    request = service.files().create(
        media_body=media,
        body={'name': file_name} # 'parents': ['<your folder Id>']}
    )
    response = None
    while response is None:
        status, response = request.next_chunk()
        if status and progress:
            print("Uploaded %d%%." % int(status.progress() * 100))
    return True

def getCredentials():
    creds = None
    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    creds_dir = os.path.join(os.path.expanduser('~'), '.gcp')
    if not os.path.exists(creds_dir):
        os.makedirs(creds_dir)
    tokenfile = os.path.join(creds_dir, 'token.pickle')
    if os.path.exists(tokenfile):
        with open(tokenfile, 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                os.path.join(creds_dir, 'credentials.json'), SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open(tokenfile, 'wb') as token:
            pickle.dump(creds, token)
    return creds

if __name__ == '__main__':
    main()
