node { // Python 2.7
	 stage('Checkout') {
        checkout scm
    }
    stage('InstallBuildTest') {
        sh 'docker build --build-arg PYTHON_VERSION=2.7 .'
    }
}

node { // Python 3.5
    stage('Checkout') {
        checkout scm
    }
    stage('InstallBuildTest') {
        sh 'docker build --build-arg PYTHON_VERSION=3.5 .'
    }
}
