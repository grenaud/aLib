{-# LANGUAGE BangPatterns, OverloadedStrings #-}
{-# LANGUAGE NamedFieldPuns, RecordWildCards #-}
{-# LANGUAGE QuasiQuotes #-}

import Control.Applicative
import Control.Monad
import Data.Binary.Strict.Get
import Data.Char
import Data.IntMap ( IntMap )
import Data.List
import Data.Maybe
import Data.Monoid
import Data.Word
import GHC.Float
import Numeric
import System.Console.GetOpt
import System.Directory
import System.Environment
import System.Exit
import System.FilePath
import System.IO
import System.Process
import Text.Blaze.Html
import Text.Blaze.Html.Renderer.Utf8
import Text.Blaze.Html4.Strict          hiding ( map, div )
import Text.Blaze.Html4.Strict.Attributes
import Text.XML.Light

import qualified Data.ByteString as S
import qualified Data.IntMap as IM
import qualified Data.Text as T
import qualified Data.Text.IO as T

data Conf = Conf { output_folder :: FilePath -> FilePath
                 , progress :: String -> IO ()
                 , debug    :: String -> IO () }

options :: [ OptDescr (Conf -> IO Conf) ]
options = [
    Option "o" ["output", "outfolder"] (ReqArg set_output "DIR")
           "Generate Report in folder DIR (def RUNFOLDER/Data/Intensities/BaseCalls/RTAreport/)",
    Option "v" ["verbose"] (NoArg set_verbose) 
           "Print more progress reports",
    Option "q" ["quiet"] (NoArg set_quiet) 
           "Print no progress reports",
    Option "h?" ["help", "usage"] (NoArg usage) 
           "Print this message"
          ]
  where
    usage _ = do pn <- getProgName
                 let blah = "Usage: " ++ pn ++ " [option...] runfolder"
                 hPutStrLn stderr $ usageInfo blah options
                 exitSuccess

    set_output fp c = return $ c { output_folder = const fp }
    set_verbose   c = return $ c { progress = hPutStrLn stderr, debug = hPutStrLn stderr }
    set_quiet     c = return $ c { progress = \_ -> return (),  debug = \_ -> return ()  }

defaultOptions :: IO Conf
defaultOptions = return $ Conf outf (hPutStrLn stderr) (\_ -> return ())
  where
    outf rf = rf </> "Data" </> "Intensities" </> "BaseCalls" </> "RTAreport"


index_length :: Int
index_length = 10

data RunInfo = RunInfo {
    run_id :: String,
    run_instrument :: String,
    run_cycles :: [Int],        -- length of reads, in order
    run_lanes :: Int,
    run_tiles :: Int }
  deriving Show

{- data RunStatus = RunStatus {
    run_version :: String
    run_flags :: M.Map String Bool }
  deriving Show -}

ranges :: [Int] -> [(Int,Int)]
ranges ls = zip (scanl (+) 1 ls) (tail $ scanl (+) 0 ls)

-- Urgh, this is very different for GA, MiSeq, HiSeq :(
-- We try to parse whatever version is thrown at us.
read_runinfo :: FilePath -> IO RunInfo
read_runinfo fp = S.readFile fp >>= \raw -> case parse raw of
    (ri:_) -> return ri
    [    ] -> fail $ "could not parse runinfo: " ++ show fp
  where 
    parse raw = do runinfo <- maybeToList $ parseXMLDoc raw
                   node <- findChildren (unqual "Run") runinfo
                   runid <- maybeToList $ findAttr (unqual "Id") node
                   instrument <- strContent <$> findChildren (unqual "Instrument") node

                   let get_length_from_range n = do l <- maybeToList (findAttr (unqual "FirstCycle") n) >>= readM
                                                    r <- maybeToList (findAttr (unqual  "LastCycle") n) >>= readM
                                                    return $ r-l+1

                   cycles <- msum [ do return $ findChildren (unqual "Reads") node >>= findChildren (unqual "Read") >>= 
                                                maybeToList . findAttr (unqual "NumCycles") >>= readM
                                  , do return $ findChildren (unqual "Reads") node >>= findChildren (unqual "Read") >>= 
                                                get_length_from_range
                                  ]
                   (lanes, tiles) <- msum [ do let ls = findChildren (unqual "Tiles") node >>= findChildren (unqual "Lane")
                                               ts <- case ls >>= maybeToList . findAttr (unqual "Incorporation") >>= readM of
                                                        [] -> mzero ; is -> return $ maximum is
                                               return (length ls, ts)

                                          , do node2 <- findChildren (unqual "FlowcellLayout") node
                                               ls <- maybeToList $ findAttr (unqual "LaneCount") node2 >>= readM
                                               ts <- maybeToList $ findAttr (unqual "TileCount") node2 >>= readM
                                               let swathes = fromMaybe 1 $ findAttr (unqual "SwathCount") node2 >>= readM
                                                   surfaces = fromMaybe 1 $ findAttr (unqual "SurfaceCount") node2 >>= readM
                                               return (ls, ts*swathes*surfaces)
                                          ]
                   return $ RunInfo runid instrument cycles lanes tiles

      {-
       -- fixup for three-read runs that actually have a second index.
       -- will skip it for now.
      res['fix']=False
      if "_PEDI" in res['runid'].upper() and len(res['read_ranges']) == 3:
        lindex = res['read_ranges'][1][1]-res['read_ranges'][1][0]+1
        res['read_ranges'][-1]=(res['read_ranges'][-1][0],res['read_ranges'][-1][1]-lindex)
        res['read_ranges'].append((res['read_ranges'][-1][1]+1,res['read_ranges'][-1][1]+lindex))
        res['fix'] = True
      -}

readM :: (Read a, Monad m, MonadPlus m) => String -> m a
readM = msum . fmap (return . fst) . filter (all isSpace . snd) . reads

-- read_status: extract RTA version number and boolean flags for
-- configuration settings {CopyAllFiles, CopyImages, DeleteImages,
-- DeleteIntensity, IsPairedEndRun} from "Status.xml".  This information
-- doesn't seem to be available any more.

-- read_stats_read: extract some values from a summary file.  Unfortunately,
-- this summary doesn't seem to exist anymore.  We may be able to
-- generate it from other sources, though.
--
-- Wanted per lane:
-- - ClustersRaw: computable from TileStatistics
-- - ClustersRawSD: computable from above
-- - ClustersPF: computable from TileStatistics
-- - ClustersPFSD: computable from above
-- - PrcPFClusters: computable from above
-- - PrcPFClustersSD: computable from above
-- - Phasing: no idea
-- - Prephasing: no idea
-- - PrcAlign: to control?  that's computable from ControlMetricsOut.bin
-- - PrcAlignSD: to control?  that's computable from above
-- - FirstCycleIntPF: huh?
-- - FirstCycleIntPFSD: huh?
-- - PrcIntensityAfter20CyclesPF: huh?
-- - PrcIntensityAfter20CyclesPFSD: huh?

{-
def det_mean_cycle_error(reads,lanes,tiles,root):
  res = []
  for (start,end) in reads:
    if (end-start+1) > index_length:
      res.append({"Start":start,'End':end,"Errors":{},"ErrorsSD":{},"AveError":{},"AveErrorSD":{}})
      for lane in range(lanes):
        res[-1]["Errors"][lane+1]=[]
        res[-1]["ErrorsSD"][lane+1]=[]
        res[-1]["AveError"][lane+1]=None
        res[-1]["AveErrorSD"][lane+1]=None
      errsum = [0]*lanes
      sdsum = [0]*lanes
      counts = [0]*lanes
      for cycle in range(start,end+1):
        # helper =^= lane ~> tile ~> TL_value
        helper = read_error_chart_xml(root+"%d.xml"%cycle)
        if helper != None:
          for lane in range(lanes):
            val = mean(helper,lane+1)
            val2 = sd(helper,lane+1,val)
            if val != None:
              errsum[lane] += val
              sdsum[lane] += val2
              counts[lane] += 1
            res[-1]["Errors"][lane+1].append(val)
            res[-1]["ErrorsSD"][lane+1].append(val2)
        else:
          for lane in range(lanes): res[-1]["Errors"][lane+1].append(None)
      for lane in range(lanes):
        res[-1]["AveError"][lane+1] = None if counts[lane] == 0 else errsum[lane]/counts[lane]
        res[-1]["AveErrorSD"][lane+1] = None if counts[lane] == 0 else sdsum[lane]/counts[lane]
  return res
-}

{-
mean_cycle_error :: (Int,Int) -> ErrorStats -> ...
mean_cycle_error ranges stats = 
    [ EStats { es_start = fromc
             , es_end   = toc
             , es_err = err
             , es_sd_err = sd_err
             , es_mean_err   = mean_err
             , es_mean_err_sd = mean_err_sd }
    | (fromc,toc) <- ranges
    , (lane, lane_stats) <- IM.toList stats
    , let (_, lane_stats1) = IM.split (fromc-1) lane_stats
    , let (lane_stats2, _) = IM.split (toc+1) lane_stats1
    , let means = map (mean . IM.elems) $ IM.elems lane_stats2
    , let sds   = zipWith sd means $ IM.elems lane_stats2
    , let count = length $ filter (not . IM.null) $ IM.elems lane_stats2
    , let ave_error = sum means `pdiv` 
    , let 

             , 

  where
    -- for a given range
    -- errsum ~ sum over cycles in range of mean error rate over tiles
    errsum = sum [ mean [ stats ! ln ! tl ! cy | ... ] | cy <- ... ]

mean_error_for_cycle :: IntMap Double -> (Double, Double
    -}
    
-- Lane ~> Cycle ~> Tile ~> ErrorRate
type ErrorStats = IntMap (IntMap (IntMap Double))

read_error_stats :: FilePath -> IO ErrorStats
read_error_stats fp = S.readFile fp >>= either fail return . fst . runGet hdr
  where
    hdr = do 3 <- getWord8
             l <- getWord8
             loop (fromIntegral l) IM.empty

    loop l !m = do e <- isEmpty
                   if e then return m
                        else do !ln <- fromIntegral <$> getWord16le
                                !tl <- fromIntegral <$> getWord16le
                                !cy <- fromIntegral <$> getWord16le
                                !v  <- float2Double <$> getFloat32host
                                skip (l-10)

                                let !m1  = IM.findWithDefault IM.empty ln m
                                    !m2  = IM.findWithDefault IM.empty cy m1
                                    !m2' = IM.insert tl v m2
                                    !m1' = IM.insert cy m2' m1
                                    !m'  = IM.insert ln m1' m
                                loop l m'


-- We get four intensity and focus values from "ExtractionMetrics".
read_intensities :: FilePath -> IO TileStatistics
read_intensities fp = S.readFile fp >>= either fail return . fst . runGet hdr
  where
    getFloat = float2Double <$> getFloat32host
    hdr = do 2 <- getWord8
             38 <- getWord8
             loop IM.empty

    add k1 k2 k3 v m = let !m2 = IM.findWithDefault IM.empty k1 m
                           !tl = IM.findWithDefault mempty k2 m2
                           !m3' = IM.insert k3 v $ tile_ext_metrics tl
                           !tl' = tl { tile_ext_metrics = m3' }
                           !m2' = IM.insert k2 tl' m2
                       in IM.insert k1 m2' m
    
    loop !ts = do e <- isEmpty
                  if e then return ts
                       else do ln <- fromIntegral <$> getWord16le
                               tl <- fromIntegral <$> getWord16le
                               cy <- fromIntegral <$> getWord16le
                               vals <- ExtMetrics <$> getFloat    <*> getFloat    <*> getFloat    <*> getFloat
                                                  <*> getWord16le <*> getWord16le <*> getWord16le <*> getWord16le
                               skip 8
                               loop $ add ln tl cy vals ts


data ExtMetrics = ExtMetrics { fwhm_a, fwhm_c, fwhm_g, fwhm_t :: !Double
                             , int_a, int_c, int_g, int_t :: !Word16 }

select_focus, select_ints :: [ExtMetrics -> Double]
select_focus = [ fwhm_a, fwhm_c, fwhm_g, fwhm_t ]
select_ints  = map ((.) fromIntegral) [ int_a, int_c, int_g, int_t ]

instance Monoid Tile where
    mempty = Tile Nothing Nothing IM.empty
    Tile a b c `mappend` Tile u v w = Tile (a `orElse` u) (b `orElse` v) (c `IM.union` w)
      where
        Nothing `orElse` b = b
        a       `orElse` _ = a


-- used to read cluster and cluster pf values; both are now in the same
-- file (tile metrics)
read_cluster_counts :: FilePath -> IO TileStatistics
read_cluster_counts fp = S.readFile fp >>= either fail return . fst . runGet hdr
  where
    hdr = do 2 <- getWord8
             10 <- getWord8
             loop IM.empty

    add k1 k2 v m = let !mm = IM.findWithDefault IM.empty k1 m
                        !v0 = IM.findWithDefault mempty k2 mm
                        !v' = v `mappend` v0
                        !m' = IM.insert k2 v' mm
                    in IM.insert k1 m' m    
    
    loop !ts = do e <- isEmpty
                  if e then return ts
                       else do ln <- fromIntegral <$> getWord16le
                               tl <- fromIntegral <$> getWord16le
                               mc <- fromIntegral <$> getWord16le
                               v  <- round <$> getFloat32host
                               case mc of
                                      102 -> loop $ add ln tl (mempty { tile_clusters    = Just v }) ts
                                      103 -> loop $ add ln tl (mempty { tile_pf_clusters = Just v }) ts
                                      _   -> loop ts


-- Lane ~> Tile ~> (density, Cycle ~> focus+intensity)
type TileStatistics = IntMap (IntMap Tile)
data Tile = Tile { tile_clusters :: Maybe Int
                 , tile_pf_clusters :: Maybe Int
                 , tile_ext_metrics :: IntMap ExtMetrics }

read_stats_html :: TileStatistics -> (Int,Int) -> String -> Html
read_stats_html stats (fromc,toc) info = do
    let headers  = [ "Lane", "Clusters [Tile]", "Clusters PF", "PF %", "1st Cycle Int (PF)"
                   , "% Intensity 20cycles", "Phasing", "Prephasing" ]

    h3 $ toMarkup $ info ++ " Read"
    table ! border "1" ! align "center" $ do
        tr $ mapM_ th headers
        sequence_ [ tr $ do td (toMarkup (show lane)) 
                            two_td 0 (fmap fromIntegral . tile_clusters) lstats
                            two_td 0 (fmap fromIntegral . tile_pf_clusters) lstats
                            two_td 2 pf_perc lstats

                            -- Since this must be "PF only", I believe
                            -- we need to get this from the "called
                            -- clusters" statistics.  Even if slightly
                            -- incorrect, it's close enough.
                            -- two_td 0 ... -- 1st cycle intensity?
                            -- two_td 2 ... -- 20th cycle perc intensity?

                            -- Phasing? Pre-Phasing?  seems to be missing.
                            -- td ! align "right" $ toMarkup $ showFFloat (Just 4) ... []
                  | (lane, lstats) <- IM.toList stats ]


  where
    pf_perc :: Tile -> Maybe Double
    pf_perc x = do p <- tile_pf_clusters x 
                   d <- tile_clusters x
                   guard $ d /= 0
                   return $ 100 * fromIntegral p / fromIntegral d

    -- Output mean and sd.  Each tile is a sample, we sum over a lane.
    two_td :: Int -> (Tile -> Maybe Double) -> IntMap Tile -> Html
    two_td prec sel xs = do
        let len  = IM.foldl' (\a x -> case sel x of Just  _ -> a + 1 ; Nothing -> a) 0 xs
            mean = IM.foldl' (\a x -> case sel x of Just xx -> a + xx ; Nothing -> a) 0 xs / len
            var  = IM.foldl' (\a x -> case sel x of Just xx -> a + (xx - mean) ^ 2 ; Nothing -> a) 0 xs / len

        td ! align "right" $ preEscapedToMarkup $
            showFFloat (Just prec) mean " &plusmn;" ++
            showFFloat (Just    2) (sqrt var) []
                                        

error_stats_html :: ErrorStats -> (Int,Int) -> String -> Html
error_stats_html errors (fromc,toc) info = do
    let middle = (toc-fromc-1) `div` 2
    let headers = [ "% Align (PF)", "Error 1st base",
                    "Error " ++ shows (middle+1) "th base",
                    "Last base", "Average Error" ]

    -- we should be able to glean this stuff from ControlMetrics.bin

    h3 $ toMarkup $ "Error of controls for " ++ info ++ " Read (" ++ shows fromc "-" ++ shows toc ")"
    p $ table ! border "1" ! align "center" $ do
        tr $ mapM_ (th . toMarkup) $ "Lane" : headers

{-
def error_stats_html(errors,read_stats,info=''):
    rrange = (errors['Start'],errors['End'])

    lanes = filter(lambda x: type(x) == type(1),read_stats.keys())
    for lane in range(min(lanes),max(lanes)+1):
      outstr += "<tr><td>%d</td>"%lane

      for cols in headers:
        if len(cols) == 3:
          if read_stats[lane][cols[1]] > 100:
            outstr += "<td align='right'>%.0f &plusmn;%.2f</td>"%(read_stats[lane][cols[1]],read_stats[lane][cols[2]])
          else:
            outstr += "<td align='right'>%.2f &plusmn;%.2f</td>"%(read_stats[lane][cols[1]],read_stats[lane][cols[2]])
        elif len(cols) == 2:
          outstr += "<td align='right'>%.4f</td>"%(read_stats[lane][cols[1]])

      for i in [0,middle,-1]:
        value = errors['Errors'][lane][i]
        if value != None: outstr += "<td align='right'>%.2f%% &plusmn;%.2f</td>"%(value,errors['ErrorsSD'][lane][i])
        else: outstr += "<td align='center'>NA</td>"

      if errors['AveError'][lane] != None: outstr += "<td align='right'>%.2f%% &plusmn;%.2f</td>"%(errors['AveError'][lane],errors['AveErrorSD'][lane])
      else: "<td align='center'>NA</td>"
-}


mean_of_im :: IntMap Double -> Double
mean_of_im m = sum (IM.elems m) / fromIntegral (IM.size m)

-- Lane ~> Cycle ~> mean error rate
type MeanErrorStats = IntMap (IntMap Double)

generate_image_error :: FilePath -> [(Int,Int)] -> MeanErrorStats -> IO Html
generate_image_error outfolder ranges errors = do
    let img_file = "error_lanes.png"

    let x_vals = concat [ [f..t] | (f,t) <- ranges ]

    let all_lanes = IM.keys errors
        maxl = shows (maximum all_lanes)

    let rcmd = unlines $
            ( "x <- c( " ++ show_list x_vals ++ ")" ) :
            ( "y <- c()" ) :

            [ "y <- cbind(y,c(" ++ show_listM errs ++ "))"
              | lane <- all_lanes
              , let m = IM.findWithDefault IM.empty lane errors
              , let errs = [ IM.lookup cy m | cy <- x_vals ] ] ++

            [ "png('" ++ (outfolder </> "images" </> img_file) ++ "',width=800,height=400)"
            , "matplot(x,y,xlab='Cycle',main='Per cycle error rate of control reads'," ++
              "ylab='Error rate [%]',type='l',lty=1,lwd=2,col=1:" ++ maxl ")"
            , "legend('topleft',sprintf('Lane %d',1:" ++ maxl "),fill=1:" ++ maxl ")"
            , "dev.off()" ]

       
    (Just hin, Nothing, Nothing, r_pid) <- withFile "/dev/null" WriteMode $ \null_hdl ->
                                           createProcess (proc "R" ["--vanilla", "--quiet"])
                                           { std_in = CreatePipe, std_out = UseHandle null_hdl }
    hPutStrLn hin rcmd
    hClose hin

    ecode <- waitForProcess r_pid
    case ecode of
        ExitSuccess   -> return $ p ! align "center" $ an_img img_file
        ExitFailure e -> fail $ "R failed with " ++ show e


an_img :: String -> Html
an_img img_file = a ! href v $ img ! src v
  where v = toValue $ "images" </> img_file

showM :: Show a => Maybe a -> String
showM Nothing  = "N/A"
showM (Just a) = show a

show_list :: Show a => [a] -> String
show_list = intercalate "," . map show

show_listM :: Show a => [Maybe a] -> String
show_listM = intercalate "," . map (maybe "NA" show)

float_or_na :: ([Double] -> Double) -> [Double] -> String
float_or_na _ [] = "N/A"
float_or_na f xs = showFFloat (Just 2) (f xs) []

pdiv :: Double -> Int -> Maybe Double
pdiv _ 0 = Nothing
pdiv a b = Just $ a / fromIntegral b

ave_ints :: ExtMetrics -> Int
ave_ints em = (sum [ fromIntegral (f em) | f <- [int_a,int_c,int_g,int_t] ] +2) `div` 4

fwhms :: ExtMetrics -> [Double]
fwhms em = [ fwhm_a em, fwhm_c em, fwhm_g em, fwhm_t em ]

show_ratio :: (Eq a, Integral a) => Maybe a -> Maybe a -> String
show_ratio (_      ) (Just  0) = "NaN"
show_ratio (Just pf) (Just cl) = showFFloat (Just 2) (100 * fromIntegral pf / fromIntegral cl :: Double) []
show_ratio _         _         = "N/A"

cluster_int_focus :: TileStatistics -> Html
cluster_int_focus ts = do
    forM_ (IM.toList ts) $ \(lane, ts') -> do
        h3 $ toMarkup $ "Lane " ++ show lane
        table ! border "1" ! align "center" $ do
            tr $ mapM_ th [ "Lane", "Tile", "Clusters", "Clusters PF", "% PF"
                          , "1st Cycle Int", "20th Cycle Int", "% Intensity"
                          , "Min Focus", "Median Focus", "Max Focus" ]
            let nrtiles = IM.size ts'
            let to_td x = td ! align "right" $ toMarkup x

            forM_ (IM.toList ts') $ \(tile, stats) ->
                tr $ do td $ toMarkup $ show lane
                        td $ toMarkup $ show tile
                        to_td $ showM $ tile_clusters stats
                        to_td $ showM $ tile_pf_clusters stats
                        to_td $ show_ratio (tile_pf_clusters stats) (tile_clusters stats)

                        let first_int  = ave_ints <$> IM.lookup  1 (tile_ext_metrics stats)
                            twenty_int = ave_ints <$> IM.lookup 21 (tile_ext_metrics stats)
                            focuses    = concatMap fwhms $ IM.elems $ tile_ext_metrics stats
                            
                        to_td $ showM first_int
                        to_td $ showM twenty_int
                        to_td $ show_ratio twenty_int first_int

                        to_td $ float_or_na minimum focuses
                        to_td $ float_or_na maximum focuses
                        to_td $ float_or_na median  focuses

median :: Ord a => [a] -> a
median xs = m xs (length xs `div` 2)
  where
    m (x:xs) k | length l == k = x
               | length l <  k = m r (k - length l - 1)
               | otherwise     = m l k
      where (l,r) = partition (<x) xs


generate_images_int_focus :: FilePath -> String -> [ExtMetrics -> Double] -> TileStatistics -> IO Html
generate_images_int_focus outfolder ctype selectors stats = do
    let lanes = IM.keys stats

    let max_cycle = maximum $ map (IM.size . tile_ext_metrics) $ concatMap IM.elems $ IM.elems stats

    let img_name lane = ctype ++ "_lane" ++ shows lane ".png"
        ctype' = case ctype of [] -> [] ; (u:v) -> toUpper u : v

    let rcmd = unlines $
            [ "x <- c(" ++ show_list [1..max_cycle] ++ ")" ] ++ concat
            [ [ base : " <- c(" ++ show_listM clist ++ ")"
              | (base, sel) <- zip "ACGT" selectors
              , let stats_lane = IM.findWithDefault IM.empty lane stats
                    per_tile_and_cycle_stat = map (IM.map sel . tile_ext_metrics) $ IM.elems stats_lane
                    ntiles = length per_tile_and_cycle_stat
                    per_cycle = IM.unionsWith (+) per_tile_and_cycle_stat
                    clist = [ IM.lookup cy per_cycle >>= (`pdiv` ntiles) | cy <- [1..max_cycle] ] 
              ] ++
              [ "png('" ++ (outfolder </> "images" </> img_name lane) ++ "',width=800,height=400)"
              , "matplot(x,cbind(A,C,G,T),xlab='Cycle',main='Per cycle " ++ ctype' ++ " values " ++
                "(Lane " ++ show lane ++ ")',ylab='" ++ ctype' ++ "',type='l',lty=1,lwd=2," ++
                "col=c('green','blue','black','red'))"
              , "legend('topright',c('A','C','G','T'),fill=c('green','blue','black','red'))"
              , "dev.off()" ]
            | lane <- lanes ]  

    (Just hin, Nothing, Nothing, r_pid) <- withFile "/dev/null" WriteMode $ \null_hdl ->
                                           createProcess (proc "R" ["--vanilla", "--quiet"])
                                           { std_in = CreatePipe, std_out = UseHandle null_hdl }
    hPutStrLn hin rcmd
    hClose hin

    ecode <- waitForProcess r_pid
    case ecode of
        ExitFailure e -> fail $ "R failed with " ++ show e
        ExitSuccess   -> return $ do h2 $ toMarkup ctype'
                                     forM_ lanes $ \lane -> do
                                         h3 $ toMarkup $ "Lane " ++ show lane
                                         p ! align "center" $ an_img $ img_name lane


main :: IO ()
main = do
    ( acts, files, errors ) <- getOpt Permute options <$> getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldl (>>=) defaultOptions acts
    runfolder <- case files of [f] -> doesDirectoryExist f >>= \e -> if e then return f else
                                        hPutStrLn stderr "runfolder does not exist" >> exitFailure
                               [ ] -> hPutStrLn stderr "missing runfolder" >> exitFailure
                               _   -> hPutStrLn stderr "exactly one runfolder expected" >> exitFailure
    let outfolder = output_folder $ runfolder
        imgfolder = outfolder </> "images"
        repfolder = runfolder </> "InterOp"

    progress "Evaluating RunInfo file..."
    RunInfo{..} <- read_runinfo $ runfolder </> "RunInfo.xml"

    -- progress "Evaluating Status file..."
    -- RunStatus{..} <- read_status $ repfolder </> "Status.xml" -- XXX

    debug "Checking for output folder..."
    createDirectoryIfMissing False outfolder
    createDirectoryIfMissing False imgfolder

    -- progress "Evaluating Read Summary files..."
    -- read_stats_forward <- read_stats_read (repfolder </> "Summary/read1.xml")

    {-
    Unfortunately, these files do no longer exists, and we'd need to generate them instead of simply copying them.

    progress "Copying some images from RTA Summary..."
    forM_ (zip ["NumClusters_Chart.png", "NumPassedFilter25_Chart.png", "PassedFilter25_Chart.png", "NumClusters By Lane.png"]
               ["num_cl_fc.png", "num_cl_pf_fc.png", "frac_cl_pf_fc.png", "pf_stats.png"]) $ \(old, new) ->
        copyFile (repfolder </> old) (outfolder </> "images" </> new)
    -}
    let pf_image = p ! align "center" $ an_img "pf_stats.png"
        flowcell_overview_images = p ! align "center" $ mapM_ an_img
                                [ "num_cl_fc.png", "num_cl_pf_fc.png", "frac_cl_pf_fc.png" ]


    progress "Determining average error rates..."
    errors <- read_error_stats (repfolder </> "ErrorMetricsOut.bin")
    -- det_mean_cycle_error(run_info['read_ranges'],run_info['lanes'],run_info['tiles'],args[0]+"/Data/reports/ErrorRate/Chart_")
    error_image <- generate_image_error outfolder (ranges run_cycles) (fmap (fmap mean_of_im) errors)

    progress "Reading tile intensity and focus values..."
    tile_stats <- read_intensities $ repfolder </> "ExtractionMetricsOut.bin" 

    intensities_dev_images <- generate_images_int_focus outfolder "intensity" select_ints tile_stats
    focus_dev_images       <- generate_images_int_focus outfolder    "focus" select_focus tile_stats

    progress "Reading tile cluster and pass filter values..."
    cl_counts <- read_cluster_counts $ repfolder </> "TileMetricsOut.bin"
 
    let eff_tile_stats = IM.unionWith (IM.unionWith mappend) tile_stats cl_counts

    progress "Generating output HTML..."

    let business_reads = filter (\(a,b) -> b-a > index_length) $ ranges run_cycles

    let report opt = docTypeHtml $ body $ do
            h1 ! align "center" $ toMarkup $ {- run_version ++ -} "Report for " ++ run_id
            p ! align "center" $ toMarkup $ intercalate ", "
                [ "Lanes: " ++ show run_lanes
                , "Tiles: " ++ show run_tiles
                , "Total number of cycles: " ++ show (sum run_cycles) ]
                
            h2 "Run Information"
            -- unfortunately not available:
            -- p $ toMarkup $ intercalate ", " [ key ++ ": " ++ if b then "True" else "False" | (key,b) <- run_flags ]
            
            h3 "Reads"
            p $ sequence_ $ intersperse br
                        [ toMarkup $ "Read " ++ show i ++ info ++ show a ++ " - " ++ show b
                        | (i, (a,b)) <- zip [1..] (ranges run_cycles)
                        , let info = if b-a <= index_length then " [Index]: " else ": " ]
 
            -- opt pf_image                        -- missing :(
            zipWithM_ (read_stats_html eff_tile_stats) business_reads ["Forward", "Reverse"]
            -- opt flowcell_overview_images        -- missing :(

            h3 "Sequencing Error"
            error_image
            zipWithM_ (error_stats_html errors) business_reads ["Forward", "Reverse"]
            
            opt $ do intensities_dev_images
                     focus_dev_images
                     h2 "Tile Statistics"
                     cluster_int_focus eff_tile_stats

    withFile (outfolder </> "Summary_short.html") WriteMode $ \hdl ->
        renderHtmlToByteStringIO (S.hPut hdl) $ report (\_ -> return ())

    withFile (outfolder </> "Summary.html") WriteMode $ \hdl ->
        renderHtmlToByteStringIO (S.hPut hdl) $ report (\x -> x)


