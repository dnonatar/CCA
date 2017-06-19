from flask import render_template, request, g, send_from_directory, Flask
from werkzeug import secure_filename
#from .forms import LoginForm
import os,sqlite3,string,random,csv,requests,glob,smtplib,re

app = Flask(__name__)

DATABASE = "/home/ratanond/Desktop/Masters_Project/CCA/Tool/tool.db" 

def connect_to_database():
    return sqlite3.connect(DATABASE)

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = connect_to_database()
        db.text_factory = str
    return db

@app.route('/')
def uploadcca():
    return render_template("upload_demonstr.html")


@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None:
        db.close()

def save_submission(query, args):
	c = get_db().cursor();
	c.execute(query, args)
	get_db().commit()
	id = c.lastrowid
	c.close()
	return id

@app.route('/importcca', methods= ['POST'])
def import_cca():
	get_db().execute("CREATE TABLE IF NOT EXISTS ccajobsss(status TEXT, name TEXT, email TEXT, project TEXT, description TEXT,rand TEXT, id integer primary key autoincrement)")
	input1 = request.files['file1']	    
	input2 = request.files['file2']
	
	#microarray_filenames = request.files['microarray_filenames']    
	#name = request.form.get('name')
	#email = request.form.get('email')
	#project = request.form.get('project')
	#description = request.form.get('description')

	if input1 and allowed_file(input1.filename):
	   dir1 = "/home/ratanond/Desktop/Masters_Project/CCA/Tool/Input1/"
	   filename = secure_filename(input1.filename)
	   input1.save(os.path.join(dir1, filename))
	    

	if file and allowed_file(input2.filename):
	   dir2 = "/home/ratanond/Desktop/Masters_Project/CCA/Tool/Input2/"
	   filename = secure_filename(input2.filename)
	   input2.save(os.path.join(dir2, filename))

	    

	
	lid = save_submission("INSERT INTO ccajobsss (status, name, email, project, description, rand) VALUES (?,?,?,?,?,?)", ["Incomplete", "name", "email", "project", "description", "abc"]);
	
	

	return render_template('job_success.html')

app.config['ALLOWED_EXTENSIONS'] = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif','csv','TXT'])

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

if __name__ == "__main__":
    app.run(debug=True)



