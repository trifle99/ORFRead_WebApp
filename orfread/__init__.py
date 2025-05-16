from flask import Flask
from .database import Base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import os

# sqlite engine for lightweight database
engine = create_engine('sqlite:///orfdatabase.db', echo=True)
Base.metadata.create_all(engine)  # create our database tabels
Session = sessionmaker(bind=engine)  # create our session

#initialize our flask app
def create_app():
    #app set up
    app=Flask(__name__)

    app.secret_key=os.environ['API_KEY']
    app.config['SECRET KEY']=os.environ['API_KEY']

    #configuring our database
    app.config['SQLALCHEMY_DATABASE_URI']='sqlite:///orfdatabase.db'

    #routes set up
    from .routes import routes
    app.register_blueprint(routes, url_prefix='/')

    return app
