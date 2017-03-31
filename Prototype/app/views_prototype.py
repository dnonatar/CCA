from flask import render_template, request, g, send_from_directory, Flask
from werkzeug import secure_filename
#from .forms import LoginForm
import os,sqlite3,string,random,csv,requests,glob,smtplib,re

app = Flask(__name__)

DATABASE = "/home/ratanond/prototype.db" 

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
    return render_template("uploadcca.html")


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
	microbiota = request.files.getlist("file[]")	    
	microarray = request.files.getlist("file[]")
	immunology = request.files.getlist("file[]")
	#microarray_filenames = request.files['microarray_filenames']    
	#name = request.form.get('name')
	#email = request.form.get('email')
	#project = request.form.get('project')
	#description = request.form.get('description')

	for file in microbiota:
	    if file and allowed_file(file.filename):
		dir = "/home/ratanond/Desktop/Thesis/Prototype/app/microbiota/"
		filename = secure_filename(file.filename)
		file.save(os.path.join(dir, filename))
	    

	for file in microarray:
	    if file and allowed_file(file.filename):
		dir = "/home/ratanond/Desktop/Thesis/Prototype/app/microarray/"
		filename = secure_filename(file.filename)
		file.save(os.path.join(dir, filename))
	    

	for file in immunology:
	    if file and allowed_file(file.filename):
		dir = "/home/ratanond/Desktop/Thesis/Prototype/app/immunology/"
		filename = secure_filename(file.filename)
		file.save(os.path.join(dir, filename))
	    

	
	lid = save_submission("INSERT INTO ccajobsss (status, name, email, project, description, rand) VALUES (?,?,?,?,?,?)", ["Incomplete", "name", "email", "project", "description", "abc"]);
	
	

	return render_template('job_success.html')

app.config['ALLOWED_EXTENSIONS'] = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif','csv','TXT'])

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

if __name__ == "__main__":
    app.run(debug=True)



