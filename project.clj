(defproject sen "1.0.0-SNAPSHOT"
  :description "Sampled ensemble neutrality"

  :jvm-opts ["-Xmx2g"]
  :main sen.robustness.neutrality
  :aot [sen.robustness.neutrality]

  :java-source-paths ["src-java/"]

  ;:offline? true
  
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [mlabs.jars/clojure-contrib "1.2.0-mlab"]
                 [org.clojure/tools.namespace "0.2.2"]
                                  
                 [slingshot "0.10.3"]   ; enhanced exceptions
                 ;;[mlabs.jars/clojure-csv "1.3.0-mlab"]
                 [clojure-csv/clojure-csv "2.0.1"]
                 [clj-shell "0.1.0"]

                 ;;command line options
                 [org.clojure/tools.cli "0.2.2"]

                 ;;[org.clojure.contrib/macro-utils "1.3.0-alpha4"]
                 [org.clojure/tools.macro "0.1.1"]
                 
                 [swingrepl "1.3.0"]
                 [incanter "1.4.0"]
                 
                 ;;file reading in parallel
                 [iota "1.1.1"]
                 [foldable-seq "0.2"]
  
                 ;; SVM stuff
                 [nz.ac.waikato.cms.weka/weka-stable "3.6.8"]
                 [tw.edu.ntu.csie/libsvm "3.1"]

                 ;;json
                 [org.clojure/data.json "0.2.1"]

                 ]
  :profiles { :dev {:dependencies [[org.apache.hadoop/hadoop-core "1.1.2"]]}}
  
  :dev-dependencies [[jline "0.9.94"]
                     ;[mlabs.jars/swank-clojure "1.5.0-sd-mlab-col-hack"]
                     ]

  :repositories
  {"commons-releases" "http://bizdirusa.com/mirrors/apache"
   "c3p0-releases"    "http://mirrors.ibiblio.org/pub/mirrors/maven2"})
