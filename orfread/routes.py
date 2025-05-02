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
        file_name=uploaded_file.filename
        uploaded_file.save(curr_dir+'\\database_entry.txt') #saves into root project dir as database_entry.txt, it will overwrite existing files, but il delete them later anyway after processed into database

        path_to_file=curr_dir+'\\database_entry.txt'

        #ReadORF methods
        from .readorf import ReadORF
        instance=ReadORF(path_to_file)

        #checking if this sequence has already been input into db
        header=instance.header()

        from . import Session
        from .database import raw_file, ORF
        session=Session()


        query=session.query(raw_file).filter_by(header=header).all()


        if query:
            flash('This file has already been submitted before, please try again with a new file.', category='error')
        else:
            instance.add_file(file_name)
            flash('File has been added to database.', category='success')

            #ADD ORF PROTEINS
            #ORF1
            instance.ORF1(curr_dir+'\\ORF1aa.txt')
            res_dict=instance.ORFSearch(curr_dir+'\\ORF1aa.txt', curr_dir+'\\ORF1prots.txt', 'F1')
            instance.add_orf(curr_dir+'\\ORF1prots.txt', res_dict, header, 1, 'forward') #for ORF1
            instance.ntd_recall(curr_dir+'\\ORF1ntds.txt', res_dict, 'F1')
            instance.add_ntd(curr_dir+'\\ORF1ntds.txt', res_dict, header)
            #ORF2
            instance.ORF2(curr_dir+'\\ORF2aa.txt')
            res_dict=instance.ORFSearch(curr_dir+'\\ORF2aa.txt', curr_dir+'\\ORF2prots.txt', 'F2')
            instance.add_orf(curr_dir+'\\ORF2prots.txt', res_dict, header, 2, 'forward') #for ORF2
            instance.ntd_recall(curr_dir+'\\ORF2ntds.txt', res_dict, 'F2')
            instance.add_ntd(curr_dir+'\\ORF2ntds.txt', res_dict, header)
            #ORF3
            instance.ORF3(curr_dir+'\\ORF3aa.txt')
            res_dict=instance.ORFSearch(curr_dir+'\\ORF3aa.txt', curr_dir+'\\ORF3prots.txt', 'F3')
            instance.add_orf(curr_dir+'\\ORF3prots.txt', res_dict, header, 3, 'forward') #for ORF3
            instance.ntd_recall(curr_dir+'\\ORF3ntds.txt', res_dict, 'F3')
            instance.add_ntd(curr_dir+'\\ORF3ntds.txt', res_dict, header)
            #mORF1
            instance.mORF1(curr_dir+'\\mORF1aa.txt')
            res_dict=instance.mORFSearch(curr_dir+'\\mORF1aa.txt', curr_dir+'\\mORF1prots.txt', 'R1')
            instance.add_orf(curr_dir+'\\mORF1prots.txt', res_dict, header, 1, 'reverse') #for mORF1
            instance.mntd_recall(curr_dir+'\\mORF1ntds.txt', res_dict, 'R1')
            instance.add_ntd(curr_dir+'\\mORF1ntds.txt', res_dict, header)
            #mORF2
            instance.mORF2(curr_dir+'\\mORF2aa.txt')
            res_dict=instance.mORFSearch(curr_dir+'\\mORF2aa.txt', curr_dir+'\\mORF2prots.txt', 'R2')
            instance.add_orf(curr_dir+'\\mORF2prots.txt', res_dict, header, 2, 'reverse') #for mORF2
            instance.mntd_recall(curr_dir+'\\mORF2ntds.txt', res_dict, 'R2')
            instance.add_ntd(curr_dir+'\\mORF2ntds.txt', res_dict, header)
            #mORF3
            instance.mORF3(curr_dir+'\\mORF3aa.txt')
            res_dict=instance.mORFSearch(curr_dir+'\\mORF3aa.txt', curr_dir+'\\mORF3prots.txt', 'R3')
            instance.add_orf(curr_dir+'\\mORF3prots.txt', res_dict, header, 3, 'reverse') #for mORF3
            instance.mntd_recall(curr_dir+'\\mORF3ntds.txt', res_dict, 'R3')
            instance.add_ntd(curr_dir+'\\mORF3ntds.txt', res_dict, header)

            #update raw_file.orf column to include how many orf proteins that file contains #could not do this before as i had to run all readorf() instance methods to get total orf proteins
            instance.update_raw_file(header)
    else:
        flash('Please submit a correct file', category='error')

    #redirects once if POST:
    return render_template('upload.html')

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
    return render_template('test.html', headers=headers, data=res, orf_headers=orf_headers, orf_data=orf_res)

if __name__=='__main__':
    print('Direct run on routes.py script')