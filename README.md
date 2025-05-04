<!-- About -->
<h1> ABOUT </h1> 
ORF Read Web App:  <br>
Web application to read raw genome sequence files. <br>
Upload your file into the web app, and it will process your file to return all ORF proteins in your sequencing file. <br>
This is saved into a database which you can query from the web app! <br>

## How to use: 
Once the web app is being run: <br>
+Go to 'UPLOAD' page and submit your sequence file for upload <br>
+Go to 'TEST' page to view your ORF protein sequences <br>

<!-- notes -->
## Note:
+This web app is still in early stage development <br>
+Intended for large whole genome sequence files. Cannot support files which contain more than one sequence YET. (eg: >sequence 1, ACTGCTCTGCTACT..., >sequence 2, ACTGCTACACA...).  
Your file should contain one sequence only (eg: >sequence 1 @ID etc, ACTGCACGATCG...). <br>

## TODO:
+Speed up database commits. Can take a while to upload large files to database...
+Add search feature to query and search for specific sequences.
