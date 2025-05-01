from flask import Flask, Blueprint, render_template, request, redirect, url_for

routes=Blueprint('routes', __name__)

##################################################
##################################################
##################################################

#home route:
@routes.route('/')
def home():

    return render_template('home.html')

##################################################

#/upload route with GET method func to render template
@routes.route('/upload', methods=['GET'])
def uploadpage():

    return render_template('upload.html')

##################################################

#/upload route with POST method to submit file data in body request
@routes.route('/upload', methods=['POST'])
def upload():
    # func to process uploaded file
    uploaded_file=request.files['file']
    if uploaded_file.filename!='': #if uploaded file name is not empty:

        import os
        curr_dir=os.getcwd()
        uploaded_file.save(curr_dir+'\\database_entry.txt') #saves into root project dir as database_entry.txt, it will overwrite existing files, but il delete them later anyway after processed into database

    #redirects once if POST:
    return redirect(url_for('routes.home'))

##################################################

#search route
@routes.route('/search', methods=['GET', 'POST'])
def search():
    if request.method=='GET':
        return render_template('search.html')

    elif request.method=='POST':
        submit_orf=request.form['ORF_PROTEIN']

        if submit_orf=='ADD':
            from .logic import DBTest
            test=DBTest('whatever')
            test.add_entry()
            return redirect(url_for('routes.test'))

        elif submit_orf=='DELETE':
            from .logic import DBTest
            test=DBTest('whatever')
            test.delete_entry()
            return redirect(url_for('routes.test'))

        else:
            return redirect(url_for('routes.home'))



##################################################

@routes.route('/test')
def test():
    from .database import raw_file
    from . import Session

    session=Session()
    res=session.query(raw_file).all()
    #return(str(res))
    headers=raw_file.__table__.columns.keys()
    return render_template('test.html', headers=headers, data=res)