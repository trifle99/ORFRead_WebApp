from flask import Flask, Blueprint, render_template, request, redirect, url_for, flash

routes=Blueprint('routes', __name__)

##################################################
##################################################
##################################################

#home route:
@routes.route('/')
def home():

    return render_template('home.html')

##################################################
##################################################
##################################################

#/upload route with GET method func to render template
@routes.route('/upload', methods=['GET'])
def uploadpage():

    return render_template('upload.html')

##################################################

#/upload route with POST method to submit file data in body request
@routes.route('/upload', methods=['POST'])
def upload():
    # func to process uploaded file>process the data with readorf.py logic to add them into our database
    uploaded_file=request.files['file']
    if uploaded_file.filename!='': #if uploaded file name is not empty:

        import os
        import glob
        curr_dir=os.getcwd()
        file_name=uploaded_file.filename
        uploaded_file.save(curr_dir+'\\database_entry.txt') #saves into root project dir as database_entry.txt, it will overwrite existing files, but il delete them later anyway after processed into database

        path_to_file=curr_dir+'\\database_entry.txt'

        #ReadORF methods
        from .readorf import ReadORF, DBcrud
        instance=ReadORF(path_to_file)
        dbcrud=DBcrud(path_to_file)
        #checking if this sequence has already been input into db
        header=instance.header()

        from . import Session
        from .database import raw_file, ORF
        session=Session()


        query=session.query(raw_file).filter_by(header=header).all()


        if query:
            flash('This file has already been submitted before, please try again with a new file.', category='error')
        else:
            dbcrud.add_file(file_name)

            flash('File has been added to database.', category='success')

            #ORF1
            instance.ORF1(curr_dir+'\\ORF1.txt')
            res=instance.newORFSearch(curr_dir+'\\ORF1.txt', curr_dir+'\\ORF1prot.csv', 'F1')
            instance.newntd_recall(curr_dir+'\\ORF1ntd.csv', res, 'F1')
            instance.combinecsv(curr_dir+'\\ORF1prot.csv', curr_dir+'\\ORF1ntd.csv', curr_dir+'\\ORF1.csv')
            if dbcrud.len_csv(curr_dir+'\\ORF1.csv')<10:
                dbcrud.simple_add(curr_dir+'\\ORF1.csv')
            elif dbcrud.len_csv(curr_dir+'\\ORF1.csv')>9:
                dbcrud.bulk_add(curr_dir+'\\ORF1.csv')

            #deleting files and cleaning up dir
            for file in glob.glob(curr_dir+'\\ORF1*.*'):
                os.remove(file)

            #ORF2
            instance.ORF2(curr_dir+'\\ORF2.txt')
            res=instance.newORFSearch(curr_dir+'\\ORF2.txt', curr_dir+'\\ORF2prot.csv', 'F2')
            instance.newntd_recall(curr_dir+'\\ORF2ntd.csv', res, 'F2')
            instance.combinecsv(curr_dir+'\\ORF2prot.csv', curr_dir+'\\ORF2ntd.csv', curr_dir+'\\ORF2.csv')
            if dbcrud.len_csv(curr_dir+'\\ORF2.csv')<10:
                dbcrud.simple_add(curr_dir+'\\ORF2.csv')
            elif dbcrud.len_csv(curr_dir+'\\ORF2.csv')>9:
                dbcrud.bulk_add(curr_dir+'\\ORF2.csv')

            #deleting files and cleaning up dir
            for file in glob.glob(curr_dir+'\\ORF2*.*'):
                os.remove(file)

            #ORF3
            instance.ORF3(curr_dir+'\\ORF3.txt')
            res=instance.newORFSearch(curr_dir+'\\ORF3.txt', curr_dir+'\\ORF3prot.csv', 'F3')
            instance.newntd_recall(curr_dir+'\\ORF3ntd.csv', res, 'F3')
            instance.combinecsv(curr_dir+'\\ORF3prot.csv', curr_dir+'\\ORF3ntd.csv', curr_dir+'\\ORF3.csv')
            if dbcrud.len_csv(curr_dir+'\\ORF3.csv')<10:
                dbcrud.simple_add(curr_dir+'\\ORF3.csv')
            elif dbcrud.len_csv(curr_dir+'\\ORF3.csv')>9:
                dbcrud.bulk_add(curr_dir+'\\ORF3.csv')

            #deleting files and cleaning up dir
            for file in glob.glob(curr_dir+'\\ORF3*.*'):
                os.remove(file)

            #nORF1
            instance.mORF1(curr_dir+'\\mORF1.txt')
            res=instance.mnewORFSearch(curr_dir+'\\mORF1.txt', curr_dir+'\\mORF1prot.csv', 'R1')
            instance.mnewntd_recall(curr_dir+'\\mORF1ntd.csv', res, 'R1')
            instance.combinecsv(curr_dir+'\\mORF1prot.csv', curr_dir+'\\mORF1ntd.csv', curr_dir+'\\mORF1.csv')
            if dbcrud.len_csv(curr_dir+'\\mORF1.csv')<10:
                dbcrud.simple_add(curr_dir+'\\mORF1.csv')
            elif dbcrud.len_csv(curr_dir+'\\mORF1.csv')>9:
                dbcrud.bulk_add(curr_dir+'\\mORF1.csv')

            #deleting files and cleaning up dir
            for file in glob.glob(curr_dir+'\\mORF1*.*'):
                os.remove(file)

            #mORF2
            instance.mORF2(curr_dir+'\\mORF2.txt')
            res=instance.mnewORFSearch(curr_dir+'\\mORF2.txt', curr_dir+'\\mORF2prot.csv', 'R2')
            instance.mnewntd_recall(curr_dir+'\\mORF2ntd.csv', res, 'R2')
            instance.combinecsv(curr_dir+'\\mORF2prot.csv', curr_dir+'\\mORF2ntd.csv', curr_dir+'\\mORF2.csv')
            if dbcrud.len_csv(curr_dir+'\\mORF2.csv')<10:
                dbcrud.simple_add(curr_dir+'\\mORF2.csv')
            elif dbcrud.len_csv(curr_dir+'\\mORF2.csv')>9:
                dbcrud.bulk_add(curr_dir+'\\mORF2.csv')

            #deleting files and cleaning up dir
            for file in glob.glob(curr_dir+'\\mORF2*.*'):
                os.remove(file)

            #mORF3
            instance.mORF3(curr_dir+'\\mORF3.txt')
            res=instance.mnewORFSearch(curr_dir+'\\mORF3.txt', curr_dir+'\\mORF3prot.csv', 'R3')
            instance.mnewntd_recall(curr_dir+'\\mORF3ntd.csv', res, 'R3')
            instance.combinecsv(curr_dir+'\\mORF3prot.csv', curr_dir+'\\mORF3ntd.csv', curr_dir+'\\mORF3.csv')
            if dbcrud.len_csv(curr_dir+'\\mORF3.csv')<10:
                dbcrud.simple_add(curr_dir+'\\mORF3.csv')
            elif dbcrud.len_csv(curr_dir+'\\mORF3.csv')>9:
                dbcrud.bulk_add(curr_dir+'\\mORF3.csv')

            #deleting files and cleaning up dir
            for file in glob.glob(curr_dir+'\\mORF3*.*'):
                os.remove(file)


            #update raw_file to show how many proteins that file has
            dbcrud.update_raw_file(header)

    else:
        flash('Please submit a correct file', category='error')

    #redirects once if POST:
    return render_template('upload.html')

##################################################
##################################################
##################################################

#search route
@routes.route('/search', methods=['GET', 'POST'])
def search():
    if request.method=='GET':
        return render_template('search.html')

    elif request.method=='POST':
        submit_orf=request.form['ORF_PROTEIN']

        if submit_orf=='ADD':
            from .readorf import DBTest
            test=DBTest('whatever')
            test.add_entry()
            return redirect(url_for('routes.test'))

        elif submit_orf=='DELETE':
            from .readorf import DBTest
            test=DBTest('whatever')
            test.delete_entry()
            return redirect(url_for('routes.test'))

        else:
            return redirect(url_for('routes.home'))

##################################################
##################################################
##################################################

@routes.route('/test')
def test():
    from .database import raw_file, ORF
    from . import Session

    session=Session()

    res=session.query(raw_file).all()
    #return(str(res))
    headers=raw_file.__table__.columns.keys()

    #ORF Table data
    orf_res=session.query(ORF).all()
    orf_headers=ORF.__table__.columns.keys()

    #paginating our ORF table data

    total=session.query(ORF).count()
    page=request.args.get('page', 1, type=int)
    per_page=50
    offset=(page-1)*per_page
    items=session.query(ORF).offset(offset).limit(per_page).all()

    pages=(total+per_page-1)//per_page

    return render_template('test.html', headers=headers, data=res, orf_headers=orf_headers, orf_data=items, page=page, pages=pages)

##################################################
##################################################
##################################################

@routes.route('/query', methods=['GET','POST'])
def query():

    if request.method=='GET':

        from .database import raw_file, ORF
        #this is to load the headers available from the database into our select form inputs
        from . import Session
        session = Session()

        from sqlalchemy import select
        # executing sql expression diretly using select and use scalars to return a scalar result object (normal result objects are returned as: [(1,),(2,)] but with scalar results they returend as [1,2]
        # distinct returns unique values within that select query
        unq_headers = session.scalars(select(ORF.ORF_header).distinct()).all()

        return render_template('query.html', select_headers=unq_headers)

    if request.method=='POST':

        #save and send form data to redirected url:query_results
        query=request.form.getlist('orf_query')
        orf_headers=request.form.getlist('orf_table')

        return redirect(url_for('routes.query_results', query_orf=query, query_headers=orf_headers))

##################################################
##################################################
##################################################

@routes.route('/query_results', methods=['GET'])
def query_results():

    if request.method=='GET':

        query_orf = request.args.getlist('query_orf')
        query_headers = request.args.getlist('query_headers')

        from .database import ORF
        from . import Session
        from sqlalchemy import select
        session_sql=Session()

        from sqlalchemy import select
        # executing sql expression diretly using select and use scalars to return a scalar result object (normal result objects are returned as: [(1,),(2,)] but with scalar results they returend as [1,2]
        # distinct returns unique values within that select query
        unq_headers = session_sql.scalars(select(ORF.ORF_header).distinct()).all()
        select_query=int(query_orf[0])

        #setting default values on form query if theyre empty
        if query_orf[1]=="":
            query_orf_len=0
        else:
            query_orf_len=query_orf[1]

        if query_orf[2]=="":
            query_ntd_len=0
        else:
            query_ntd_len=query_orf[2]

        if query_orf[3]=="":
            query_ntd_ind1=0
        else:
            query_ntd_ind1=query_orf[3]

        if query_orf[4]=='':
            #query for last ntd_index

            query_ntd_ind_last=session_sql.query(ORF).order_by(-ORF.NTD_index).limit(1).first()
            import re
            #re to search for ntd_index number
            pattern=r'NTD_index:(\d+)'
            match=re.search(pattern, str(query_ntd_ind_last))
            query_ntd_ind2=int(match.group(1))+1

        else:
            query_ntd_ind2=query_orf[4]

        #add pagination to query data
        page=request.args.get('page', 1, type=int) #getting value of page arg in url, default value if no page in url is 1 with type integer
        per_page=50 #how many items we want to display per page
        offset=(page-1)*per_page #starting point/index of each page
        total=session_sql.query(ORF).filter(ORF.ORF_header==unq_headers[select_query], ORF.ORF_len>query_orf_len, ORF.NTD_len>query_ntd_len, ORF.NTD_index>query_ntd_ind1, ORF.NTD_index<query_ntd_ind2).count() #total count of items in our filtered data
        pages=(total+per_page-1)//per_page #total number of pages

        #adding query form filters to our query and adding offset+limit for paginations:
        query_data=session_sql.query(ORF).filter(ORF.ORF_header==unq_headers[select_query], ORF.ORF_len>query_orf_len, ORF.NTD_len>query_ntd_len, ORF.NTD_index>query_ntd_ind1, ORF.NTD_index<query_ntd_ind2).offset(offset).limit(per_page).all() #our data


    #i have to pass our form_data (query_headers, query_orf) in to html <a href...> during pagination so that it remembers the form_data filters
    return render_template('query_results.html', query_headers=query_headers, query_orf=query_orf, query_data=query_data, page=page, pages=pages)

##################################################
##################################################
##################################################

if __name__=='__main__':

    print('Direct run on routes.py script')