<!-- About -->
<h1> ABOUT </h1> 
ORF Read Web App:  <br>
Web application to read raw genome sequence files. <br>
Upload your file into the web app, and it will process your file to return all ORF proteins in your sequencing file. <br>
This is saved into a database which you can query from the web app! <br>

## How to use: 
Once the web app is being run: <br>
+Go to 'UPLOAD' page and submit your sequence file for upload <br>
+Go to 'QUERY' page and select which ORF file you would like to view your ORF sequences for <br>
+If you want to search for certain features, you can apply filters by editing the QUERY forms <br>
+If you want to view the whole database, just leave the QUERY form blank <br>

<!-- notes -->
## Note:
+This web app is still in early stage development <br>
+Intended for large whole genome sequence files. Cannot support files which contain more than one sequence YET. (eg: >sequence 1, ACTGCTCTGCTACT..., >sequence 2, ACTGCTACACA...).  
Your file should contain one sequence only (eg: >sequence 1 @ID etc, ACTGCACGATCG...). <br>

## TODO:
+Speed up database commits. Can take a while to upload large files to database... :white_check_mark: before: took ~30mins to add 40,000 entries into ORF, now it takes ~10s to add 140,000 entries <br>
+Add search feature to query and search for specific sequences. <br> :white_check_mark: can now query database with certain filters by filling out query form <br>

