node { // Python 2.7
	 stage('Checkout') {
        checkout scm
    }
    stage('PythonEnvironmentSetup') {
        sh 'wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -O miniconda.sh'
        sh 'bash miniconda.sh -b -p .'
    }
    stage('Install') {
        sh 'PATH=$PWD/miniconda/bin:$PATH bash devtools/ci/jenkins/install.sh'
    }
    stage('Test') {
        sh 'PATH=$PWD/miniconda/bin:$PATH bash devtools/ci/jenkins/runtest.sh'
    }
    stage('Cleanup') {
        sh 'bash devtools/ci/jenkins/cleanup.sh'
    }
}
