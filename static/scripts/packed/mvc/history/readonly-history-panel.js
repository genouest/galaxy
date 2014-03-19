define(["mvc/history/history-model","mvc/dataset/hda-base"],function(h,a){var g=SessionStorageModel.extend({defaults:{expandedHdas:{},show_deleted:false,show_hidden:false},addExpandedHda:function(j){var i="expandedHdas";this.save(i,_.extend(this.get(i),_.object([j],[true])))},removeExpandedHda:function(j){var i="expandedHdas";this.save(i,_.omit(this.get(i),j))},toString:function(){return"HistoryPrefs("+this.id+")"}});g.storageKeyPrefix="history:";g.historyStorageKey=function f(i){if(!i){throw new Error("HistoryPrefs.historyStorageKey needs valid id: "+i)}return(g.storageKeyPrefix+i)};g.get=function d(i){return new g({id:g.historyStorageKey(i)})};g.clearAll=function c(j){for(var i in sessionStorage){if(i.indexOf(g.storageKeyPrefix)===0){sessionStorage.removeItem(i)}}};var b=Backbone.View.extend(LoggableMixin).extend({HDAViewClass:a.HDABaseView,tagName:"div",className:"history-panel",fxSpeed:"fast",emptyMsg:_l("This history is empty"),noneFoundMsg:_l("No matching datasets found"),initialize:function(i){i=i||{};if(i.logger){this.logger=i.logger}this.log(this+".initialize:",i);this.linkTarget=i.linkTarget||"_blank";this.fxSpeed=_.has(i,"fxSpeed")?(i.fxSpeed):(this.fxSpeed);this.filters=[];this.searchFor="";this.findContainerFn=i.findContainerFn;this.hdaViews={};this.indicator=new LoadingIndicator(this.$el);this._setUpListeners();var j=_.pick(i,"initiallyExpanded","show_deleted","show_hidden");this.setModel(this.model,j,false);if(i.onready){i.onready.call(this)}},_setUpListeners:function(){this.on("error",function(j,m,i,l,k){this.errorHandler(j,m,i,l,k)});this.on("loading-history",function(){this._showLoadingIndicator("loading history...",40)});this.on("loading-done",function(){this._hideLoadingIndicator(40);if(_.isEmpty(this.hdaViews)){this.trigger("empty-history",this)}});this.once("rendered",function(){this.trigger("rendered:initial",this);return false});if(this.logger){this.on("all",function(i){this.log(this+"",arguments)},this)}return this},errorHandler:function(k,n,j,m,l){console.error(k,n,j,m,l);if(n&&n.status===0&&n.readyState===0){}else{if(n&&n.status===502){}else{var i=this._parseErrorMessage(k,n,j,m,l);if(!this.$messages().is(":visible")){this.once("rendered",function(){this.displayMessage("error",i.message,i.details)})}else{this.displayMessage("error",i.message,i.details)}}}},_parseErrorMessage:function(l,p,k,o,n){var j=Galaxy.currUser,i={message:this._bePolite(o),details:{user:(j instanceof User)?(j.toJSON()):(j+""),source:(l instanceof Backbone.Model)?(l.toJSON()):(l+""),xhr:p,options:(p)?(_.omit(k,"xhr")):(k)}};_.extend(i.details,n||{});if(p&&_.isFunction(p.getAllResponseHeaders)){var m=p.getAllResponseHeaders();m=_.compact(m.split("\n"));m=_.map(m,function(q){return q.split(": ")});i.details.xhr.responseHeaders=_.object(m)}return i},_bePolite:function(i){i=i||_l("An error occurred while getting updates from the server");return i+". "+_l("Please contact a Galaxy administrator if the problem persists.")},loadHistoryWithHDADetails:function(k,j,i,m){var l=function(n){return _.keys(g.get(n.id).get("expandedHdas"))};return this.loadHistory(k,j,i,m,l)},loadHistory:function(l,k,j,o,m){var i=this;k=k||{};i.trigger("loading-history",i);var n=h.History.getHistoryData(l,{historyFn:j,hdaFn:o,hdaDetailIds:k.initiallyExpanded||m});return i._loadHistoryFromXHR(n,k).fail(function(r,p,q){i.trigger("error",i,r,k,_l("An error was encountered while "+p),{historyId:l,history:q||{}})}).always(function(){i.trigger("loading-done",i)})},_loadHistoryFromXHR:function(k,j){var i=this;k.then(function(l,m){i.JSONToModel(l,m,j)});k.fail(function(m,l){i.render()});return k},JSONToModel:function(l,i,j){this.log("JSONToModel:",l,i,j);j=j||{};if(Galaxy&&Galaxy.currUser){l.user=Galaxy.currUser.toJSON()}var k=new h.History(l,i,j);this.setModel(k);return this},setModel:function(j,i,k){i=i||{};k=(k!==undefined)?(k):(true);this.log("setModel:",j,i,k);this.freeModel();this.selectedHdaIds=[];if(j){if(Galaxy&&Galaxy.currUser){j.user=Galaxy.currUser.toJSON()}this.model=j;if(this.logger){this.model.logger=this.logger}this._setUpWebStorage(i.initiallyExpanded,i.show_deleted,i.show_hidden);this._setUpModelEventHandlers();this.trigger("new-model",this)}if(k){this.render()}return this},freeModel:function(){if(this.model){this.model.clearUpdateTimeout();this.stopListening(this.model);this.stopListening(this.model.hdas)}this.freeHdaViews();return this},freeHdaViews:function(){this.hdaViews={};return this},_setUpWebStorage:function(j,i,k){this.storage=new g({id:g.historyStorageKey(this.model.get("id"))});if(_.isObject(j)){this.storage.set("exandedHdas",j)}if(_.isBoolean(i)){this.storage.set("show_deleted",i)}if(_.isBoolean(k)){this.storage.set("show_hidden",k)}this.trigger("new-storage",this.storage,this);this.log(this+" (init'd) storage:",this.storage.get());return this},_setUpModelEventHandlers:function(){this.model.hdas.on("add",this.addHdaView,this);this.model.on("error error:hdas",function(j,l,i,k){this.errorHandler(j,l,i,k)},this);return this},render:function(k,l){this.log("render:",k,l);k=(k===undefined)?(this.fxSpeed):(k);var i=this,j;if(this.model){j=this.renderModel()}else{j=this.renderWithoutModel()}$(i).queue("fx",[function(m){if(k&&i.$el.is(":visible")){i.$el.fadeOut(k,m)}else{m()}},function(m){i.$el.empty();if(j){i.$el.append(j.children())}m()},function(m){if(k&&!i.$el.is(":visible")){i.$el.fadeIn(k,m)}else{m()}},function(m){if(l){l.call(this)}i.trigger("rendered",this);m()}]);return this},renderWithoutModel:function(){var i=$("<div/>"),j=$("<div/>").addClass("message-container").css({"margin-left":"4px","margin-right":"4px"});return i.append(j)},renderModel:function(){var i=$("<div/>");i.append(b.templates.historyPanel(this.model.toJSON()));this.$emptyMessage(i).text(this.emptyMsg);i.find(".history-secondary-actions").prepend(this._renderSearchButton());this._setUpBehaviours(i);this.renderHdas(i);return i},_renderEmptyMsg:function(k){var j=this,i=j.$emptyMessage(k);if(!_.isEmpty(j.hdaViews)){i.hide()}else{if(j.searchFor){i.text(j.noneFoundMsg).show()}else{i.text(j.emptyMsg).show()}}return this},_renderSearchButton:function(i){return faIconButton({title:_l("Search datasets"),classes:"history-search-btn",faIcon:"fa-search"})},_setUpBehaviours:function(i){i=i||this.$el;i.find("[title]").tooltip({placement:"bottom"});this._setUpSearchInput(i.find(".history-search-controls .history-search-input"));return this},$container:function(){return(this.findContainerFn)?(this.findContainerFn.call(this)):(this.$el.parent())},$datasetsList:function(i){return(i||this.$el).find(".datasets-list")},$messages:function(i){return(i||this.$el).find(".message-container")},$emptyMessage:function(i){return(i||this.$el).find(".empty-history-message")},renderHdas:function(j){j=j||this.$el;var i=this,l={},k=this.model.hdas.getVisible(this.storage.get("show_deleted"),this.storage.get("show_hidden"),this.filters);this.$datasetsList(j).empty();if(k.length){k.each(function(n){var m=n.get("id"),o=i._createHdaView(n);l[m]=o;if(_.contains(i.selectedHdaIds,m)){o.selected=true}i.attachHdaView(o.render(),j)})}this.hdaViews=l;this._renderEmptyMsg(j);return this.hdaViews},_createHdaView:function(j){var i=j.get("id"),k=new this.HDAViewClass({model:j,linkTarget:this.linkTarget,expanded:this.storage.get("expandedHdas")[i],hasUser:this.model.ownedByCurrUser(),logger:this.logger});this._setUpHdaListeners(k);return k},_setUpHdaListeners:function(j){var i=this;j.on("error",function(l,n,k,m){i.errorHandler(l,n,k,m)});j.on("body-expanded",function(k){i.storage.addExpandedHda(k)});j.on("body-collapsed",function(k){i.storage.removeExpandedHda(k)});return this},attachHdaView:function(k,j){j=j||this.$el;var i=this.$datasetsList(j);i.prepend(k.$el);return this},addHdaView:function(l){this.log("add."+this,l);var j=this;if(!l.isVisible(this.storage.get("show_deleted"),this.storage.get("show_hidden"))){return j}$({}).queue([function k(n){var m=j.$emptyMessage();if(m.is(":visible")){m.fadeOut(j.fxSpeed,n)}else{n()}},function i(m){var n=j._createHdaView(l);j.hdaViews[l.id]=n;n.render().$el.hide();j.scrollToTop();j.attachHdaView(n);n.$el.slideDown(j.fxSpeed)}]);return j},refreshHdas:function(j,i){if(this.model){return this.model.refresh(j,i)}return $.when()},events:{"click .message-container":"clearMessages","click .history-search-btn":"toggleSearchControls"},collapseAllHdaBodies:function(){_.each(this.hdaViews,function(i){i.toggleBodyVisibility(null,false)});this.storage.set("expandedHdas",{});return this},toggleShowDeleted:function(i){i=(i!==undefined)?(i):(!this.storage.get("show_deleted"));this.storage.set("show_deleted",i);this.renderHdas();return this.storage.get("show_deleted")},toggleShowHidden:function(i){i=(i!==undefined)?(i):(!this.storage.get("show_hidden"));this.storage.set("show_hidden",i);this.renderHdas();return this.storage.get("show_hidden")},_setUpSearchInput:function(j){var k=this,l=".history-search-input";function i(m){if(k.model.hdas.haveDetails()){k.searchHdas(m);return}k.$el.find(l).searchInput("toggle-loading");k.model.hdas.fetchAllDetails({silent:true}).always(function(){k.$el.find(l).searchInput("toggle-loading")}).done(function(){k.searchHdas(m)})}j.searchInput({initialVal:k.searchFor,name:"history-search",placeholder:"search datasets",classes:"history-search",onfirstsearch:i,onsearch:_.bind(this.searchHdas,this),onclear:_.bind(this.clearHdaSearch,this)});return j},toggleSearchControls:function(k,i){var j=this.$el.find(".history-search-controls"),l=(jQuery.type(k)==="number")?(k):(this.fxSpeed);i=(i!==undefined)?(i):(!j.is(":visible"));if(i){j.slideDown(l,function(){$(this).find("input").focus()})}else{j.slideUp(l)}return i},searchHdas:function(i){var j=this;this.searchFor=i;this.filters=[function(k){return k.matchesAll(j.searchFor)}];this.trigger("search:searching",i,this);this.renderHdas();return this},clearHdaSearch:function(i){this.searchFor="";this.filters=[];this.trigger("search:clear",this);this.renderHdas();return this},_showLoadingIndicator:function(j,i,k){i=(i!==undefined)?(i):(this.fxSpeed);if(!this.indicator){this.indicator=new LoadingIndicator(this.$el,this.$el.parent())}if(!this.$el.is(":visible")){this.indicator.show(0,k)}else{this.$el.fadeOut(i);this.indicator.show(j,i,k)}},_hideLoadingIndicator:function(i,j){i=(i!==undefined)?(i):(this.fxSpeed);if(this.indicator){this.indicator.hide(i,j)}},displayMessage:function(n,o,m){var k=this;this.scrollToTop();var l=this.$messages(),i=$("<div/>").addClass(n+"message").html(o);if(!_.isEmpty(m)){var j=$('<a href="javascript:void(0)">Details</a>').click(function(){Galaxy.modal.show(k._messageToModalOptions(n,o,m));return false});i.append(" ",j)}return l.html(i)},_messageToModalOptions:function(m,o,l){var i=this,n=$("<div/>"),k={title:"Details"};function j(p){p=_.omit(p,_.functions(p));return["<table>",_.map(p,function(r,q){r=(_.isObject(r))?(j(r)):(r);return'<tr><td style="vertical-align: top; color: grey">'+q+'</td><td style="padding-left: 8px">'+r+"</td></tr>"}).join(""),"</table>"].join("")}if(_.isObject(l)){k.body=n.append(j(l))}else{k.body=n.html(l)}k.buttons={Ok:function(){Galaxy.modal.hide();i.clearMessages()}};return k},clearMessages:function(){this.$messages().empty();return this},scrollPosition:function(){return this.$container().scrollTop()},scrollTo:function(i){this.$container().scrollTop(i);return this},scrollToTop:function(){this.$container().scrollTop(0);return this},scrollToId:function(j){if((!j)||(!this.hdaViews[j])){return this}var i=this.hdaViews[j];this.scrollTo(i.el.offsetTop);return this},scrollToHid:function(i){var j=this.model.hdas.getByHid(i);if(!j){return this}return this.scrollToId(j.id)},toString:function(){return"ReadOnlyHistoryPanel("+((this.model)?(this.model.get("name")):(""))+")"}});var e=['<div class="history-controls">','<div class="history-search-controls">','<div class="history-search-input"></div>',"</div>",'<div class="history-title">',"<% if( history.name ){ %>",'<div class="history-name"><%= history.name %></div>',"<% } %>","</div>",'<div class="history-subtitle clear">',"<% if( history.nice_size ){ %>",'<div class="history-size"><%= history.nice_size %></div>',"<% } %>",'<div class="history-secondary-actions"></div>',"</div>","<% if( history.deleted ){ %>",'<div class="warningmessagesmall"><strong>',_l("You are currently viewing a deleted history!"),"</strong></div>","<% } %>",'<div class="message-container">',"<% if( history.message ){ %>",'<div class="<%= history.status %>message"><%= history.message %></div>',"<% } %>","</div>",'<div class="quota-message errormessage">',_l("You are over your disk quota."),_l("Tool execution is on hold until your disk usage drops below your allocated quota."),"</div>",'<div class="tags-display"></div>','<div class="annotation-display"></div>','<div class="history-dataset-actions">','<div class="btn-group">','<button class="history-select-all-datasets-btn btn btn-default"','data-mode="select">',_l("All"),"</button>",'<button class="history-deselect-all-datasets-btn btn btn-default"','data-mode="select">',_l("None"),"</button>","</div>",'<button class="history-dataset-action-popup-btn btn btn-default">',_l("For all selected"),"...</button>","</div>","</div>",'<div class="datasets-list"></div>','<div class="empty-history-message infomessagesmall">',_l("Your history is empty. Click 'Get Data' on the left pane to start"),"</div>"].join("");b.templates={historyPanel:function(i){return _.template(e,i,{variable:"history"})}};return{ReadOnlyHistoryPanel:b}});